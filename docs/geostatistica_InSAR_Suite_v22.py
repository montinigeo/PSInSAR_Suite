"""
geostatistica_InSAR_Suite_v22.py
================================
Script geostatistica estratto da InSAR Suite v2.2.
Conservato per sviluppo futuro del modulo geostatistica.

Contiene:
  - Calcolo semivariogramma sperimentale isotropo e direzionale
  - Fit modelli teorici (sferico, esponenziale, gaussiano)
  - Analisi anisotropia con superficie polare ed ellisse
  - Finestra interattiva a sei schede con pannello parametri a due colonne
    (Max cont. / Min cont.) con nugget, sill e range per ciascuna direzione
  - Kriging ordinario con anisotropia geometrica (pykrige)
  - Cross-validation leave-one-out
  - Verifica di normalità delle velocità PS (scheda integrata)
  - Scelta del campo da interpolare (velocità TS o campo numerico del layer)

Utilizzo dalla console Python di QGIS:
    exec(open(r'C:\percorso\kriging_interattivo.py').read())

Oppure dalla console standalone (senza QGIS):
    python kriging_interattivo.py --test

In modalità standalone genera dati sintetici per il test.

Funzionalità aggiuntive rispetto a geostatistica.py v2.1:
  - Finestra di parametri MODIFICABILE dopo il calcolo automatico
  - Aggiustamento manuale di: modello, nugget, sill, range
  - Aggiustamento manuale di: direzione e rapporto di anisotropia
  - Ricalcolo kriging + aggiornamento grafici e CV on-the-fly
  - Controllo minimo punti (>= 10) con messaggio chiaro
"""

import sys
import re
import tempfile
import numpy as np
import pandas as pd

# ── Rileva ambiente (QGIS o standalone) ─────────────────────────────────────
from qgis.PyQt.QtWidgets import (
    QApplication, QDialog, QVBoxLayout, QHBoxLayout, QWidget,
    QTabWidget, QLabel, QPushButton, QComboBox, QDoubleSpinBox,
    QSpinBox, QGroupBox, QDialogButtonBox, QMessageBox,
    QTextEdit, QFileDialog, QSizePolicy, QSplitter,
    QRadioButton, QFrame
)
from qgis.PyQt.QtCore import QVariant
from qgis.PyQt.QtCore import Qt
from qgis.utils import iface
from qgis.core import (QgsProject, QgsRasterLayer,
                        QgsCoordinateReferenceSystem, QgsTask,
                        QgsApplication)
from qgis.gui import QgsProjectionSelectionDialog
IN_QGIS = True

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from scipy.optimize import curve_fit
from scipy.interpolate import griddata
from scipy.spatial.distance import pdist

try:
    from pykrige.ok import OrdinaryKriging
    HAS_PYKRIGE = True
except ImportError:
    HAS_PYKRIGE = False

try:
    from osgeo import gdal, osr
    HAS_GDAL = True
except ImportError:
    HAS_GDAL = False

try:
    from pyproj import Transformer
    HAS_PYPROJ = True
except ImportError:
    HAS_PYPROJ = False



# ============================================================
# DIALOGO SCELTA CAMPO
# ============================================================
class SceltaCampoDialog(QDialog):
    """
    Prima finestra: l'utente sceglie il valore da interpolare con il kriging.
    Due opzioni:
      A) Velocità calcolata dalla regressione lineare sulle serie DYYYYMMDD
      B) Campo numerico già presente nel layer (es. VEL_LOS, E, U, pc_mov...)
    """
    def __init__(self):
        super().__init__(None)
        self.setWindowTitle('InSAR TS – Geostatistica: scegli il valore da interpolare')
        self.setMinimumWidth(420)
        layout = QVBoxLayout(self)

        lbl = QLabel('Seleziona il valore da utilizzare per variogramma e kriging:')
        lbl.setWordWrap(True)
        layout.addWidget(lbl)

        grp = QGroupBox()
        grp.setFlat(True)
        v = QVBoxLayout(grp)

        self.rb_ts = QRadioButton(
            'Velocità calcolata dalla serie storica (regressione lineare su campi DYYYYMMDD)'
        )
        self.rb_ts.setChecked(True)
        v.addWidget(self.rb_ts)

        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        v.addWidget(line)

        self.rb_field = QRadioButton('Campo numerico del layer:')
        v.addWidget(self.rb_field)

        self.cb_field = QComboBox()
        self.cb_field.setEnabled(False)
        v.addWidget(self.cb_field)

        layout.addWidget(grp)
        self.rb_ts.toggled.connect(lambda checked: self.cb_field.setEnabled(not checked))

        btn = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btn.accepted.connect(self.accept)
        btn.rejected.connect(self.reject)
        layout.addWidget(btn)

    def populate(self, layer):
        self.cb_field.clear()
        from qgis.PyQt.QtCore import QVariant
        numeric_qvariant = {
            QVariant.Int, QVariant.UInt, QVariant.LongLong, QVariant.ULongLong,
            QVariant.Double
        }
        for field in layer.fields():
            if field.type() in numeric_qvariant:
                self.cb_field.addItem(field.name())

    def getChoice(self):
        if self.rb_ts.isChecked():
            return ('ts', None)
        else:
            return ('field', self.cb_field.currentText())

# ============================================================
# SUPERFICIE POLARE NORMALIZZATA
# ============================================================
def interpola_polare(ang, lags, gamma):
    from scipy.interpolate import griddata
    ang2 = np.concatenate([ang, ang + 180])
    g2   = np.concatenate([gamma, gamma], axis=0)
    ang2 = np.mod(ang2, 360)
    a_ext = np.concatenate([ang2 - 360, ang2, ang2 + 360])
    g_ext = np.concatenate([g2, g2, g2], axis=0)
    T, R  = np.meshgrid(a_ext, lags, indexing='ij')
    Z     = g_ext
    m     = ~np.isnan(Z)
    pts   = np.column_stack((T[m], R[m]))
    vals  = Z[m]
    if len(vals) < 4:
        return None, None, None
    ti = np.linspace(0, 360, 360)
    ri = np.linspace(lags.min(), lags.max(), 200)
    Tg, Rg = np.meshgrid(ti, ri)
    Zi = griddata(pts, vals, (Tg, Rg), method="cubic")
    vmax = np.nanmax(Zi)
    if vmax and vmax > 0:
        Zi /= vmax
    return Tg, Rg, Zi


def _make_polare_figure(ang, h_dir, g_dir, i_max, i_min,
                        ang_min_manual=None, ang_max_manual=None,
                        r_max_val=None, r_min_val=None):
    """
    Crea la figura polare.
    - La superficie (Tg, Rg, Zi) è empirica e non cambia con i parametri.
    - Le linee blu/rossa usano ang_min_manual e ang_max_manual.
    - Se r_max_val e r_min_val sono forniti, disegna l'ellisse dell'anisotropia
      con i due assi e la curva ellittica in nero.
    """
    Tg, Rg, Zi = interpola_polare(ang, h_dir, g_dir)
    if Tg is None:
        return None
    # Direzioni da mostrare
    a_min = ang_min_manual if ang_min_manual is not None else ang[i_min]
    a_max = ang_max_manual if ang_max_manual is not None else ang[i_max]
    rmax = np.nanmax(Rg)

    fig = Figure(figsize=(7, 7), tight_layout=True)
    ax  = fig.add_subplot(111, projection='polar')
    ax.pcolormesh(np.radians(Tg), Rg, Zi, shading='auto', cmap='RdYlGn_r')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)

    # ── Linee direzioni ──────────────────────────────────────────────────
    ax.plot([np.radians(a_max), np.radians(a_max + 180)],
            [rmax, rmax], color='blue', lw=2,
            label=f'Max continuità: {a_max:.1f}°')
    ax.plot([np.radians(a_min), np.radians(a_min + 180)],
            [rmax, rmax], color='red', lw=2,
            label=f'Min continuità: {a_min:.1f}°')

    # ── Ellisse anisotropia ──────────────────────────────────────────────
    if r_max_val is not None and r_min_val is not None and r_min_val > 0:
        # r_max_val e r_min_val sono in metri reali, sulla stessa scala
        # dei lag del variogramma (asse radiale del grafico = rmax metri).
        # Disegniamo direttamente in unità metriche senza normalizzare.
        phi = np.radians(a_max)  # asse maggiore = dir massima continuità

        # Ellisse in coordinate polari:
        # r(t) = r_max * r_min / sqrt((r_min*cos(t-phi))^2 + (r_max*sin(t-phi))^2)
        t_arr = np.linspace(0, 2 * np.pi, 360)
        dtheta = t_arr - phi
        denom = np.sqrt((r_min_val * np.cos(dtheta))**2 +
                        (r_max_val * np.sin(dtheta))**2)
        r_ell = r_max_val * r_min_val / np.where(denom > 0, denom, 1e-9)
        # Clip al bordo del grafico per non uscire dall'asse radiale
        r_ell = np.minimum(r_ell, rmax)

        # Ellisse — linea nera continua
        ax.plot(t_arr, r_ell, color='black', lw=1.5,
                label=f'Ellisse (r_max={r_max_val:.0f} m, r_min={r_min_val:.0f} m)')

        # Asse maggiore (range max) — segmento nero tratteggiato
        ax.plot([phi, phi + np.pi], [r_max_val, r_max_val],
                color='black', lw=1.2, ls='--', alpha=0.8)

        # Asse minore (range min) — segmento nero punteggiato
        phi_min = phi + np.pi / 2
        ax.plot([phi_min, phi_min + np.pi], [r_min_val, r_min_val],
                color='black', lw=1.2, ls=':', alpha=0.8)

    ax.set_title("Superficie polare normalizzata", y=1.08, fontsize=10)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=2, fontsize=7)
    return fig

# ============================================================
# MODELLI VARIOGRAMMA
# ============================================================
MODELLI = {
    'exponential': lambda h, nug, sill, r: nug + sill * (1 - np.exp(-h / r)),
    'spherical':   lambda h, nug, sill, r: np.where(
        h <= r,
        nug + sill * (1.5 * (h / r) - 0.5 * (h / r) ** 3),
        nug + sill
    ),
    'gaussian':    lambda h, nug, sill, r: nug + sill * (1 - np.exp(-(h / r) ** 2)),
}

# Modelli disponibili: spherical, exponential, gaussian
# Modelli disponibili: spherical, exponential, gaussian

MIN_PUNTI = 10  # soglia minima punti per il kriging


# ============================================================
# VERIFICA DI NORMALITÀ
# ============================================================
def _make_normalita_figure(values, label='Valori'):
    """
    Produce una figura con 4 grafici:
    - Istogramma con curva normale teorica sovrapposta
    - Q-Q plot
    - Boxplot
    - Statistiche + risultato Shapiro-Wilk
    """
    from scipy import stats

    v = values[~np.isnan(values)]
    n = len(v)

    fig = Figure(figsize=(12, 9), tight_layout=True)
    fig.suptitle(f"Verifica di Normalità — {label}  (n={n})", fontsize=12, y=0.98)

    # ── Istogramma + curva normale ───────────────────────────────────
    ax1 = fig.add_subplot(2, 2, 1)
    mu, sigma = np.mean(v), np.std(v)
    ax1.hist(v, bins=min(20, n // 3 + 1), density=True,
             color='steelblue', edgecolor='black', alpha=0.7, label='Dati')
    xmin, xmax = ax1.get_xlim()
    x = np.linspace(xmin, xmax, 200)
    ax1.plot(x, stats.norm.pdf(x, mu, sigma), 'r-', lw=2, label='Normale teorica')
    ax1.set_title("Istogramma + Normale teorica")
    ax1.set_xlabel(label)
    ax1.set_ylabel("Densità")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3)

    # ── Q-Q plot ─────────────────────────────────────────────────────
    ax2 = fig.add_subplot(2, 2, 2)
    (osm, osr_), (slope, intercept, r) = stats.probplot(v, dist="norm")
    ax2.scatter(osm, osr_, color='steelblue', s=20, alpha=0.7, zorder=3)
    x_line = np.array([min(osm), max(osm)])
    ax2.plot(x_line, slope * x_line + intercept, 'r-', lw=2, label=f'R²={r**2:.4f}')
    ax2.set_title("Q-Q Plot (Normale)")
    ax2.set_xlabel("Quantili teorici")
    ax2.set_ylabel("Quantili osservati")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)

    # ── Boxplot ──────────────────────────────────────────────────────
    ax3 = fig.add_subplot(2, 2, 3)
    bp = ax3.boxplot(v, vert=True, patch_artist=True,
                     boxprops=dict(facecolor='steelblue', alpha=0.7),
                     medianprops=dict(color='red', lw=2))
    ax3.set_title("Boxplot")
    ax3.set_ylabel(label)
    ax3.grid(alpha=0.3, axis='y')

    # ── Statistiche + Shapiro-Wilk ───────────────────────────────────
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis('off')

    skewness = stats.skew(v)
    kurt     = stats.kurtosis(v)  # excess kurtosis

    # Shapiro-Wilk (affidabile per n < 5000)
    if n <= 5000:
        W, p_sw = stats.shapiro(v)
        sw_txt  = f"W = {W:.4f}\np-value = {p_sw:.4f}"
        if p_sw < 0.05:
            sw_esito = "⚠ Rifiuta H₀ (p<0.05)\nDistribuzione non normale"
            sw_color = '#C62828'
        else:
            sw_esito = "✓ Non rifiuta H₀ (p≥0.05)\nCompatibile con normalità"
            sw_color = '#2E7D32'
    else:
        W, p_sw  = np.nan, np.nan
        sw_txt   = "n > 5000: Shapiro-Wilk\nnon applicabile"
        sw_esito = "(usa Q-Q plot e istogramma)"
        sw_color = '#555555'

    stats_txt = (
        f"Media:          {mu:.4f}\n"
        f"Dev. Std:       {sigma:.4f}\n"
        f"Min:            {np.min(v):.4f}\n"
        f"Max:            {np.max(v):.4f}\n"
        f"Skewness:       {skewness:.4f}\n"
        f"Kurtosis (exc): {kurt:.4f}\n"
        f"N punti:        {n}\n"
        f"\n── Shapiro-Wilk ──\n"
        f"{sw_txt}"
    )
    ax4.text(0.05, 0.97, stats_txt, transform=ax4.transAxes,
             fontsize=9, verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='#F5F5F5', alpha=0.8))
    ax4.text(0.05, 0.22, sw_esito, transform=ax4.transAxes,
             fontsize=9, verticalalignment='top', color=sw_color,
             bbox=dict(boxstyle='round', facecolor='#FAFAFA', alpha=0.8))
    ax4.set_title("Statistiche")

    return fig

# ============================================================
# CALCOLO VARIOGRAMMI
# ============================================================
def _coppie(coords, v, max_d):
    MAX_PAIRS = 500_000
    n = len(v)
    n_pairs = n * (n - 1) // 2
    if n_pairs > MAX_PAIRS:
        rng = np.random.default_rng(42)
        n_sample = int(np.sqrt(2 * MAX_PAIRS)) + 1
        idx = rng.choice(n, size=min(n_sample, n), replace=False)
        coords, v = coords[idx], v[idx]

    h_all = pdist(coords, metric='euclidean')
    i_idx, j_idx = np.triu_indices(len(v), k=1)
    dv = v[i_idx] - v[j_idx]
    g_all = 0.5 * dv ** 2
    dx_all = coords[j_idx, 0] - coords[i_idx, 0]
    dy_all = coords[j_idx, 1] - coords[i_idx, 1]
    mask = h_all <= max_d
    return h_all[mask], g_all[mask], dx_all[mask], dy_all[mask]


def semivariogramma_isotropo(coords, v, n_lags, max_d):
    h, g, _, _ = _coppie(coords, v, max_d)
    bins = np.linspace(0, max_d, n_lags + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])
    gamma = np.array([
        np.mean(g[(h >= bins[k]) & (h < bins[k + 1])]) if np.any((h >= bins[k]) & (h < bins[k + 1])) else np.nan
        for k in range(n_lags)
    ])
    return centers, gamma


def semivariogrammi_direzionali(coords, v, n_lags, max_d, ang_step):
    h_all, g_all, dx_all, dy_all = _coppie(coords, v, max_d)
    ang = np.mod(np.degrees(np.arctan2(dy_all, dx_all)), 180)
    bins_ang = np.arange(0, 180 + ang_step, ang_step)
    ang_c = 0.5 * (bins_ang[:-1] + bins_ang[1:])
    bins_h = np.linspace(0, max_d, n_lags + 1)
    h_c = 0.5 * (bins_h[:-1] + bins_h[1:])
    gamma = np.full((len(ang_c), n_lags), np.nan)
    for i in range(len(ang_c)):
        m_ang = (ang >= bins_ang[i]) & (ang < bins_ang[i + 1])
        if not np.any(m_ang):
            continue
        h_dir, g_dir = h_all[m_ang], g_all[m_ang]
        for k in range(n_lags):
            m = (h_dir >= bins_h[k]) & (h_dir < h_c[k])
            if np.any(m):
                gamma[i, k] = np.mean(g_dir[m])
    return ang_c, h_c, gamma


def fit_variogram(h, g, modello=None):
    m = ~np.isnan(g)
    if np.sum(m) < 4:
        return None
    h, g = h[m], g[m]
    p0 = [np.min(g) * 0.1, np.max(g) - np.min(g), h[len(h) // 2]]
    results = {}
    candidati = {modello: MODELLI[modello]} if modello else MODELLI
    for name, f in candidati.items():
        try:
            p, _ = curve_fit(f, h, g, p0=p0,
                             bounds=([0, 0, 0.1], [np.inf, np.inf, np.max(h)]))
            e = np.mean((g - f(h, *p)) ** 2)
            results[name] = (p, e)
        except Exception:
            pass
    if not results:
        return None
    # Preferisce sferico se il suo MSE è entro il 20% del migliore
    best_name = min(results, key=lambda k: results[k][1])
    best_err  = results[best_name][1]
    if 'spherical' in results:
        sph_err = results['spherical'][1]
        if sph_err <= best_err * 1.20:
            best_name = 'spherical'
    return (best_name, results[best_name][0])


# ============================================================
# KRIGING E CROSS-VALIDATION
# ============================================================
def kriging_ordinario(coords, values, variogram_fit,
                      anisotropy_angle=None, anisotropy_ratio=None,
                      n_cells=150):
    if not HAS_PYKRIGE:
        raise ImportError("pykrige non disponibile")
    x, y = coords[:, 0], coords[:, 1]
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)
    dx, dy = xmax - xmin, ymax - ymin
    pixel_size = max(dx, dy) / (n_cells - 1)
    nx = int(np.round(dx / pixel_size)) + 1
    ny = int(np.round(dy / pixel_size)) + 1
    grid_x = np.linspace(xmin, xmin + (nx - 1) * pixel_size, nx)
    grid_y = np.linspace(ymin, ymin + (ny - 1) * pixel_size, ny)
    nugget, sill, vrange = variogram_fit[1]
    OK = OrdinaryKriging(
        x, y, values,
        variogram_model=variogram_fit[0],
        variogram_parameters={'nugget': nugget, 'sill': sill, 'range': vrange},
        anisotropy_scaling=anisotropy_ratio if anisotropy_ratio and anisotropy_ratio > 1 else 1.0,
        anisotropy_angle=np.radians(anisotropy_angle) if anisotropy_angle is not None else 0.0,
        verbose=False, enable_plotting=False
    )
    z, ss = OK.execute('grid', grid_x, grid_y)
    return grid_x, grid_y, z, ss, pixel_size


def cross_validation(coords, values, variogram_fit,
                     anisotropy_angle=None, anisotropy_ratio=None,
                     max_cv=150):
    if not HAS_PYKRIGE:
        raise ImportError("pykrige non disponibile")
    n = len(values)
    rng = np.random.default_rng(42)
    idx_cv = rng.choice(n, size=min(max_cv, n), replace=False)
    predicted, observed = np.full(n, np.nan), np.full(n, np.nan)
    nugget, sill, vrange = variogram_fit[1]
    for i in idx_cv:
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        if np.isnan(values[i]):
            continue
        observed[i] = values[i]
        try:
            OK = OrdinaryKriging(
                coords[mask, 0], coords[mask, 1], values[mask],
                variogram_model=variogram_fit[0],
                variogram_parameters={'nugget': nugget, 'sill': sill, 'range': vrange},
                anisotropy_scaling=anisotropy_ratio if anisotropy_ratio and anisotropy_ratio > 1 else 1.0,
                anisotropy_angle=np.radians(anisotropy_angle) if anisotropy_angle is not None else 0.0,
                verbose=False, enable_plotting=False
            )
            z, _ = OK.execute('points',
                              np.array([coords[i, 0]]), np.array([coords[i, 1]]))
            predicted[i] = z[0]
        except Exception:
            pass
    valid = ~np.isnan(observed) & ~np.isnan(predicted)
    res = observed[valid] - predicted[valid]
    return {
        'residuals': res,
        'observed': observed[valid],
        'predicted': predicted[valid],
        'rmse': np.sqrt(np.mean(res ** 2)),
        'mae': np.mean(np.abs(res)),
        'bias': np.mean(res),
        'std': np.std(res),
        'n_points': len(res)
    }


# ============================================================
# GRAFICI
# ============================================================
def _make_variogram_figure(h, g, fit, title):
    fig = Figure(figsize=(6, 4), tight_layout=True)
    ax = fig.add_subplot(111)
    valid = ~np.isnan(g)
    ax.scatter(h[valid], g[valid], c='k', s=30, zorder=3)
    if fit:
        name, p = fit
        hh = np.linspace(0, np.nanmax(h), 300)
        ax.plot(hh, MODELLI[name](hh, *p), 'r', lw=2,
                label=f"{name}\nnug={p[0]:.2f} sill={p[1]:.2f} r={p[2]:.0f}")
        ax.legend(fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Distanza h [m]")
    ax.set_ylabel("Semivarianza γ(h)")
    ax.grid(alpha=0.3)
    return fig


def _make_kriging_figure(grid_x, grid_y, z, coords, values):
    fig = Figure(figsize=(7, 5), tight_layout=True)
    ax = fig.add_subplot(111)
    im = ax.imshow(np.flipud(z), origin='upper',
                   extent=(grid_x.min(), grid_x.max(),
                           grid_y.min(), grid_y.max()),
                   cmap='RdYlGn', aspect='auto')
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=values,
                    cmap='RdYlGn', edgecolors='k', linewidths=0.5,
                    s=40, zorder=5)
    fig.colorbar(im, ax=ax, label="[mm/anno]")
    ax.set_title("Mappa Kriging Ordinario")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    return fig


def _make_cv_figure(cv):
    fig = Figure(figsize=(11, 4), tight_layout=True)
    ax1 = fig.add_subplot(121)
    ax1.hist(cv['residuals'], bins=15, edgecolor='black', alpha=0.7, color='skyblue')
    ax1.axvline(0, color='red', ls='--', lw=1.5)
    ax1.axvline(cv['bias'], color='orange', lw=2,
                label=f"Bias: {cv['bias']:.3f}")
    ax1.set_title("Distribuzione errori CV")
    ax1.set_xlabel("Errore [mm/anno]")
    ax1.set_ylabel("Frequenza")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3)
    ax2 = fig.add_subplot(122)
    mn = min(cv['observed'].min(), cv['predicted'].min())
    mx = max(cv['observed'].max(), cv['predicted'].max())
    ax2.scatter(cv['observed'], cv['predicted'], alpha=0.7, color='steelblue')
    ax2.plot([mn, mx], [mn, mx], 'r--', lw=2)
    r = np.corrcoef(cv['observed'], cv['predicted'])[0, 1]
    ax2.text(0.05, 0.93, f"R={r:.3f}  RMSE={cv['rmse']:.3f}",
             transform=ax2.transAxes,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7),
             fontsize=8)
    ax2.set_title("Osservato vs Predetto")
    ax2.set_xlabel("Osservato [mm/anno]")
    ax2.set_ylabel("Predetto [mm/anno]")
    ax2.grid(alpha=0.3)
    ax2.axis('equal')
    fig.suptitle("Cross-Validation", fontsize=12, y=0.98)
    return fig


# ============================================================
# FINESTRA PRINCIPALE INTERATTIVA
# ============================================================
class KrigingInterattivoDialog(QDialog):
    """
    Finestra con pannello parametri modificabili + grafici aggiornati in tempo reale.
    Layout:
      ┌─────────────────────────────────────────────────────┐
      │  [Pannello parametri]  │  [Tab grafici]             │
      │                        │                            │
      │  Modello               │  Variogramma isotropo      │
      │  Nugget / Sill / Range │  Variogramma max/min cont. │
      │  Dir. anisotropia      │  Mappa Kriging             │
      │  Rapporto anisotropia  │  Cross-Validation          │
      │                        │                            │
      │  [Ricalcola]           │  [Parametri numerici]      │
      │  [Salva GeoTIFF]       │                            │
      └─────────────────────────────────────────────────────┘
    """
    def __init__(self, coords, values, n_lags=12, ang_step=10,
                 epsg=32632, label='Velocità [mm/anno]', parent=None):
        super().__init__(parent)
        self.setWindowTitle("InSAR Suite – Kriging Interattivo v2.2")
        self.resize(1400, 800)

        self.coords  = coords
        self.values  = values
        self.n_lags  = n_lags
        self.ang_step = ang_step
        self.epsg    = epsg
        self.label   = label

        # Calcolo iniziale automatico
        self._calc_variogrammi()
        self._auto_fit()

        # Layout principale con splitter
        main = QHBoxLayout(self)
        splitter = QSplitter(Qt.Horizontal)
        main.addWidget(splitter)

        # ── Pannello sinistro: parametri ──────────────────────────────
        left = QWidget()
        lv = QVBoxLayout(left)
        lv.setSpacing(6)

        lv.addWidget(QLabel("<b>Parametri Variogramma</b>"))

        # Modello
        grp_mod = QGroupBox("Modello variogramma isotropo equivalente")
        gm = QVBoxLayout(grp_mod)
        self.cb_modello = QComboBox()
        self.cb_modello.addItems(["spherical", "exponential", "gaussian"])
        self.cb_modello.setCurrentText("spherical")
        gm.addWidget(self.cb_modello)
        lv.addWidget(grp_mod)

        # Nugget / Sill / Range
        # ── Direzione anisotropia ──────────────────────────────────────────
        grp_dir = QGroupBox("Direzione anisotropia")
        gd = QVBoxLayout(grp_dir)
        self.sp_dir   = self._spin("Dir. min. continuità (°)", 0, 180, 1)
        self.lbl_dir_max = QLabel("Dir. max. continuità: —")
        self.lbl_dir_max.setStyleSheet("color: #1565C0; padding-left: 4px;")
        gd.addWidget(self.sp_dir)
        gd.addWidget(self.lbl_dir_max)
        if self._aniso_angle is not None:
            self.sp_dir.setValue(float(self._aniso_angle))
            self._aggiorna_lbl_dir_max(float(self._aniso_angle))
        self.sp_dir.valueChanged.connect(
            lambda v: self._aggiorna_lbl_dir_max(v))
        lv.addWidget(grp_dir)

        # ── Parametri direzionali (due colonne) ────────────────────────────
        grp_params = QGroupBox("Parametri variogramma")
        gp = QVBoxLayout(grp_params)

        # Header colonne
        hdr = QHBoxLayout()
        lbl_hdr0 = QLabel("")
        lbl_hdr0.setFixedWidth(60)
        lbl_max = QLabel("<b>Max cont.</b>")
        lbl_max.setAlignment(Qt.AlignCenter)
        lbl_min = QLabel("<b>Min cont.</b>")
        lbl_min.setAlignment(Qt.AlignCenter)
        hdr.addWidget(lbl_hdr0)
        hdr.addWidget(lbl_max)
        hdr.addWidget(lbl_min)
        gp.addLayout(hdr)

        def _dspin(vmin, vmax, dec):
            sp = QDoubleSpinBox()
            sp.setRange(vmin, vmax)
            sp.setDecimals(dec)
            return sp

        # Nugget
        row_nug = QHBoxLayout()
        lbl_n = QLabel("Nugget"); lbl_n.setFixedWidth(60)
        self.sp_nugget_max = _dspin(0, 1e6, 3)
        self.sp_nugget_min = _dspin(0, 1e6, 3)
        row_nug.addWidget(lbl_n)
        row_nug.addWidget(self.sp_nugget_max)
        row_nug.addWidget(self.sp_nugget_min)
        gp.addLayout(row_nug)

        # Sill
        row_sill = QHBoxLayout()
        lbl_s = QLabel("Sill"); lbl_s.setFixedWidth(60)
        self.sp_sill_max = _dspin(0, 1e6, 3)
        self.sp_sill_min = _dspin(0, 1e6, 3)
        row_sill.addWidget(lbl_s)
        row_sill.addWidget(self.sp_sill_max)
        row_sill.addWidget(self.sp_sill_min)
        gp.addLayout(row_sill)

        # Range
        row_range = QHBoxLayout()
        lbl_r = QLabel("Range (m)"); lbl_r.setFixedWidth(60)
        self.sp_range_max = _dspin(1, 1e7, 0)
        self.sp_range_max.setMaximum(self._lags_max)
        self.sp_range_min = _dspin(1, 1e7, 0)
        self.sp_range_min.setMaximum(self._lags_max)
        row_range.addWidget(lbl_r)
        row_range.addWidget(self.sp_range_max)
        row_range.addWidget(self.sp_range_min)
        gp.addLayout(row_range)

        lv.addWidget(grp_params)

        # Nota kriging
        lbl_nota = QLabel(
            "ℹ  Kriging usa i parametri Max cont. (range max = range isotropo equivalente)"
        )
        lbl_nota.setStyleSheet("color: #888; font-size: 10px; padding: 2px;")
        lbl_nota.setWordWrap(True)
        lv.addWidget(lbl_nota)

        # Inizializza valori
        self._init_param_values()

        # Parametri variogramma sperimentale
        grp_lag = QGroupBox("Variogramma sperimentale")
        gl = QVBoxLayout(grp_lag)
        self.sp_nlags = QSpinBox()
        self.sp_nlags.setRange(4, 40)
        self.sp_nlags.setValue(n_lags)
        self.sp_nlags.setPrefix("N. lag: ")
        self.sp_angstep = QSpinBox()
        self.sp_angstep.setRange(5, 45)
        self.sp_angstep.setValue(ang_step)
        self.sp_angstep.setSingleStep(5)
        self.sp_angstep.setPrefix("Step angolare: ")
        self.sp_angstep.setSuffix("°")
        gl.addWidget(self.sp_nlags)
        gl.addWidget(self.sp_angstep)
        lv.addWidget(grp_lag)

        # Pulsanti
        btn_calc = QPushButton("🔄  Ricalcola")
        btn_calc.setStyleSheet("font-weight:bold; padding:6px;")
        btn_calc.clicked.connect(self._ricalcola)
        lv.addWidget(btn_calc)

        btn_auto = QPushButton("⚡  Ripristina automatico")
        btn_auto.clicked.connect(self._ripristina_auto)
        lv.addWidget(btn_auto)

        if HAS_GDAL:
            btn_tif = QPushButton("💾  Salva GeoTIFF...")
            btn_tif.clicked.connect(self._salva_geotiff)
            lv.addWidget(btn_tif)

        lv.addStretch()

        btn_close = QPushButton("Chiudi")
        btn_close.clicked.connect(self.accept)
        lv.addWidget(btn_close)

        left.setMaximumWidth(280)
        splitter.addWidget(left)

        # ── Pannello destro: grafici ───────────────────────────────────
        right = QWidget()
        rv = QVBoxLayout(right)
        self.tabs = QTabWidget()
        rv.addWidget(self.tabs)
        splitter.addWidget(right)
        splitter.setSizes([280, 1100])

        # Pannello testo parametri
        self.txt_params = QTextEdit()
        self.txt_params.setReadOnly(True)
        self.txt_params.setMaximumHeight(120)
        rv.addWidget(self.txt_params)

        # Prima visualizzazione
        self._aggiorna_grafici()

    # ── Helpers ──────────────────────────────────────────────────────────
    def _spin(self, prefix, vmin, vmax, decimals):
        sp = QDoubleSpinBox()
        sp.setRange(vmin, vmax)
        sp.setDecimals(decimals)
        sp.setPrefix(f"{prefix}: ")
        sp.setSingleStep(10 ** (-decimals) * 10)
        sp.setMinimumWidth(200)
        return sp

    def _calc_variogrammi(self):
        """Ricalcola variogrammi sperimentali con i parametri correnti."""
        md = np.nanmax(np.linalg.norm(
            self.coords[:, None] - self.coords[None], axis=2))
        self._md = md
        self._h_iso, self._g_iso = semivariogramma_isotropo(
            self.coords, self.values, self.n_lags, md)
        self._ang, self._h_dir, self._g_dir = semivariogrammi_direzionali(
            self.coords, self.values, self.n_lags, md, self.ang_step)
        # Distanza massima effettiva del grafico polare = centro ultimo lag
        self._lags_max = float(self._h_dir.max()) if len(self._h_dir) > 0 else md

    def _auto_fit(self):
        """Fit automatico dei parametri."""
        self._fit_iso = fit_variogram(self._h_iso, self._g_iso)

        fits = [fit_variogram(self._h_dir, self._g_dir[i])
                for i in range(len(self._ang))]
        ranges = np.array([f[1][2] if f else np.nan for f in fits])

        i_max = np.nanargmax(ranges)
        ang_ortho = (self._ang[i_max] + 90) % 180
        tol = self.ang_step / 2
        cand = [i for i, a in enumerate(self._ang)
                if abs((a - ang_ortho + 90) % 180 - 90) < tol] or list(range(len(self._ang)))
        cand_r = [(i, ranges[i]) for i in cand if not np.isnan(ranges[i])]
        i_min = (min(cand_r, key=lambda x: x[1])[0]
                 if cand_r else (i_max + len(self._ang) // 2) % len(self._ang))

        self._fits     = fits
        self._i_max    = i_max
        self._i_min    = i_min

        if fits[i_max] and fits[i_min] and ranges[i_max] > 0:
            self._aniso_angle  = self._ang[i_min]
            self._aniso_ratio  = ranges[i_max] / ranges[i_min]
            # Range direzionali reali dai variogrammi direzionali
            self._range_max_dir = float(min(ranges[i_max], self._lags_max))
            self._range_min_dir = float(min(ranges[i_min], self._lags_max))
        else:
            self._aniso_angle   = None
            self._aniso_ratio   = 1.0
            r_iso = self._fit_iso[1][2] if self._fit_iso else self._lags_max * 0.5
            self._range_max_dir = float(min(r_iso, self._lags_max))
            self._range_min_dir = float(min(r_iso, self._lags_max))

        self._last_fit = self._fit_iso

        # Pre-calcola la superficie polare (costosa, fatta una volta sola)
        self._Tg, self._Rg, self._Zi = interpola_polare(
            self._ang, self._h_dir, self._g_dir)

    def _leggi_parametri(self):
        """Legge i parametri dalla GUI.
        Il kriging usa i parametri della direzione di massima continuità
        come variogramma isotropo equivalente.
        """
        modello = self.cb_modello.currentText()
        r_max = self.sp_range_max.value()
        r_min = self.sp_range_min.value()
        # Variogramma isotropo equivalente = parametri dir. massima continuità
        p = np.array([
            self.sp_nugget_max.value(),
            self.sp_sill_max.value(),
            r_max
        ])
        aniso_angle = self.sp_dir.value()
        aniso_ratio = r_max / r_min if r_min > 0 else 1.0
        return (modello, p), aniso_angle, aniso_ratio

    def _ricalcola(self):
        """Ricalcola variogrammi se i parametri lag sono cambiati, poi aggiorna."""
        self.n_lags   = self.sp_nlags.value()
        self.ang_step = self.sp_angstep.value()
        self._calc_variogrammi()
        # Aggiorna il limite degli spinbox range con il nuovo lags_max
        self.sp_range_max.setMaximum(self._lags_max)
        self.sp_range_min.setMaximum(self._lags_max)
        # Ricalcola la superficie polare solo quando cambiano i lag
        self._Tg, self._Rg, self._Zi = interpola_polare(
            self._ang, self._h_dir, self._g_dir)
        self._aggiorna_grafici()

    def _init_param_values(self):
        """Inizializza i valori degli spinbox dai fit direzionali automatici."""
        # Max continuità: dai fit direzionali
        if self._fits and self._i_max < len(self._fits) and self._fits[self._i_max]:
            p_max = self._fits[self._i_max][1]
            self.sp_nugget_max.setValue(float(p_max[0]))
            self.sp_sill_max.setValue(float(p_max[1]))
            self.sp_range_max.setValue(float(min(p_max[2], self._lags_max)))
        elif self._fit_iso:
            p = self._fit_iso[1]
            self.sp_nugget_max.setValue(float(p[0]))
            self.sp_sill_max.setValue(float(p[1]))
            self.sp_range_max.setValue(float(min(p[2], self._lags_max)))
        # Min continuità: dai fit direzionali
        if self._fits and self._i_min < len(self._fits) and self._fits[self._i_min]:
            p_min = self._fits[self._i_min][1]
            self.sp_nugget_min.setValue(float(p_min[0]))
            self.sp_sill_min.setValue(float(p_min[1]))
            self.sp_range_min.setValue(float(min(p_min[2], self._lags_max)))
        elif self._fit_iso:
            p = self._fit_iso[1]
            self.sp_nugget_min.setValue(float(p[0]))
            self.sp_sill_min.setValue(float(p[1]))
            self.sp_range_min.setValue(float(min(p[2], self._lags_max)))

    def _aggiorna_lbl_dir_max(self, ang_min):
        """Aggiorna la label della direzione di massima continuità."""
        ang_max = (ang_min + 90) % 180
        self.lbl_dir_max.setText(f"Dir. max. continuità: {ang_max:.1f}°")

    def _ripristina_auto(self):
        """Ripristina i parametri calcolati automaticamente."""
        self._auto_fit()
        self.cb_modello.setCurrentText("spherical")
        if self._aniso_angle is not None:
            self.sp_dir.setValue(float(self._aniso_angle))
            self._aggiorna_lbl_dir_max(float(self._aniso_angle))
        self._init_param_values()
        self._aggiorna_grafici()

    def _aggiorna_grafici(self):
        """Ricalcola kriging e CV con i parametri correnti e aggiorna i tab."""
        fit, aniso_angle, aniso_ratio = self._leggi_parametri()
        self._last_fit = fit
        self._last_aniso_angle = aniso_angle
        self._last_aniso_ratio = aniso_ratio

        # Kriging
        try:
            gx, gy, z, _, ps = kriging_ordinario(
                self.coords, self.values, fit, aniso_angle, aniso_ratio)
            self._last_raster = (gx, gy, z, ps)
        except Exception as e:
            QMessageBox.warning(self, "Errore kriging", str(e))
            return

        # Cross-validation
        try:
            cv = cross_validation(self.coords, self.values, fit,
                                  aniso_angle, aniso_ratio)
        except Exception as e:
            cv = None

        # Costruisce grafici
        fig_norm   = _make_normalita_figure(self.values, self.label)
        # La superficie è pre-calcolata; solo le linee usano i parametri manuali
        ang_min_man = self.sp_dir.value()
        ang_max_man = (ang_min_man + 90) % 180
        fig_polare = _make_polare_figure(
            self._ang, self._h_dir, self._g_dir,
            self._i_max, self._i_min,
            ang_min_manual=ang_min_man,
            ang_max_manual=ang_max_man,
            r_max_val=self.sp_range_max.value(),
            r_min_val=self.sp_range_min.value())
        fig_iso = _make_variogram_figure(
            self._h_iso, self._g_iso, fit, "Variogramma Isotropo")
        # Fit direzionali con i parametri delle rispettive colonne
        modello_scelto = fit[0]
        r_max  = self.sp_range_max.value()
        r_min  = self.sp_range_min.value()
        # Grafici direzionali: usano i parametri delle rispettive colonne
        fit_max_man = (modello_scelto, np.array([
            self.sp_nugget_max.value(), self.sp_sill_max.value(), r_max]))
        fit_min_man = (modello_scelto, np.array([
            self.sp_nugget_min.value(), self.sp_sill_min.value(), r_min]))
        fig_max = _make_variogram_figure(
            self._h_dir, self._g_dir[self._i_max],
            fit_max_man,
            f"Massima Continuità ({self._ang[self._i_max]:.1f}°)")
        fig_min = _make_variogram_figure(
            self._h_dir, self._g_dir[self._i_min],
            fit_min_man,
            f"Minima Continuità ({self._ang[self._i_min]:.1f}°)")
        fig_krig = _make_kriging_figure(gx, gy, z, self.coords, self.values)
        fig_cv   = _make_cv_figure(cv) if cv else None

        # Aggiorna tab
        current = self.tabs.currentIndex()
        self.tabs.clear()
        for fig, name in [
            (fig_norm,   "Normalità"),
            (fig_polare, "Superficie Polare"),
            (fig_iso,    "Variogramma Isotropo"),
            (fig_max,    "Massima Continuità"),
            (fig_min,    "Minima Continuità"),
        ]:
            w = QWidget()
            v = QVBoxLayout(w)
            v.setContentsMargins(0, 0, 0, 0)
            canvas = FigureCanvas(fig)
            v.addWidget(canvas)
            btn_png = QPushButton("Salva PNG")
            btn_png.setFixedWidth(100)
            btn_png.clicked.connect(lambda _, f=fig, n=name: self._salva_png(f, n))
            v.addWidget(btn_png)
            self.tabs.addTab(w, name)

        # ── Tab speciale "Mappa Kriging" con pulsanti aggiuntivi ──────────
        w_krig = QWidget()
        v_krig = QVBoxLayout(w_krig)
        v_krig.setContentsMargins(0, 0, 0, 0)
        canvas_krig = FigureCanvas(fig_krig)
        v_krig.addWidget(canvas_krig)
        row_krig = QHBoxLayout()
        btn_png_krig = QPushButton("Salva PNG")
        btn_png_krig.setFixedWidth(100)
        btn_png_krig.clicked.connect(
            lambda _, f=fig_krig: self._salva_png(f, "Mappa_Kriging"))
        row_krig.addWidget(btn_png_krig)
        if IN_QGIS and HAS_GDAL and hasattr(self, '_last_raster'):
            btn_qgis = QPushButton("📥  Carica in QGIS come layer temporaneo")
            btn_qgis.setStyleSheet("font-weight: bold; padding: 4px 10px;")
            btn_qgis.clicked.connect(self._carica_in_qgis)
            row_krig.addWidget(btn_qgis)
        row_krig.addStretch()
        v_krig.addLayout(row_krig)
        self.tabs.addTab(w_krig, "Mappa Kriging")

        if fig_cv:
            w = QWidget()
            v = QVBoxLayout(w)
            v.setContentsMargins(0, 0, 0, 0)
            canvas = FigureCanvas(fig_cv)
            v.addWidget(canvas)
            btn_png_cv = QPushButton("Salva PNG")
            btn_png_cv.setFixedWidth(100)
            btn_png_cv.clicked.connect(
                lambda _, f=fig_cv: self._salva_png(f, "Cross_Validation"))
            v.addWidget(btn_png_cv)
            self.tabs.addTab(w, "Cross-Validation")

        self.tabs.setCurrentIndex(min(current, self.tabs.count() - 1))

        # Aggiorna testo parametri
        r_max_val = self.sp_range_max.value()
        r_min_val = self.sp_range_min.value()
        ratio_calc = r_max_val / r_min_val if r_min_val > 0 else 1.0
        txt = (f"Modello: {fit[0]}  |  Dir. min. continuità: {aniso_angle:.1f}°  |  "
               f"Rapporto anisotropia: {ratio_calc:.2f}\n"
               f"Max cont. — Nugget: {self.sp_nugget_max.value():.3f}  |  "
               f"Sill: {self.sp_sill_max.value():.3f}  |  "
               f"Range: {r_max_val:.0f} m  (= range isotropo equivalente)\n"
               f"Min cont. — Nugget: {self.sp_nugget_min.value():.3f}  |  "
               f"Sill: {self.sp_sill_min.value():.3f}  |  "
               f"Range: {r_min_val:.0f} m")
        if cv:
            txt += (f"\nCV — RMSE: {cv['rmse']:.4f}  |  "
                    f"MAE: {cv['mae']:.4f}  |  "
                    f"Bias: {cv['bias']:.4f}  |  "
                    f"N punti: {cv['n_points']}")
        self.txt_params.setText(txt)

    def _carica_in_qgis(self):
        """Scrive un GeoTIFF temporaneo e lo carica in QGIS come layer temporaneo."""
        if not hasattr(self, '_last_raster'):
            QMessageBox.warning(self, "Errore", "Nessun raster disponibile.")
            return
        gx, gy, z, ps = self._last_raster
        try:
            tmp = tempfile.NamedTemporaryFile(
                suffix='.tif', prefix='kriging_', delete=False)
            tmp_path = tmp.name
            tmp.close()
            _scrivi_geotiff(tmp_path, gx, gy, z, ps, self.epsg)
            layer_name = f"Kriging_{self.label} [temp]"
            rl = QgsRasterLayer(tmp_path, layer_name)
            if rl.isValid():
                QgsProject.instance().addMapLayer(rl)
                QMessageBox.information(
                    self, "Layer caricato",
                    f"Il raster è stato caricato in QGIS come:\n'{layer_name}'\n\n"
                    "È un layer temporaneo — per salvarlo permanentemente\n"
                    "usa il tasto destro sul layer → Esporta → Salva come..."
                )
            else:
                QMessageBox.warning(self, "Errore",
                    "Impossibile caricare il raster in QGIS.")
        except Exception as e:
            QMessageBox.warning(self, "Errore", str(e))

    def _salva_png(self, fig, name):
        fname, _ = QFileDialog.getSaveFileName(
            self, "Salva PNG", name.replace(" ", "_") + ".png",
            "PNG (*.png)")
        if fname:
            fig.savefig(fname, dpi=300, bbox_inches='tight', facecolor='white')

    def _salva_geotiff(self):
        if not hasattr(self, '_last_raster'):
            QMessageBox.warning(self, "Errore", "Nessun raster disponibile.")
            return
        fname, _ = QFileDialog.getSaveFileName(
            self, "Salva GeoTIFF", "Kriging.tif", "GeoTIFF (*.tif)")
        if not fname:
            return
        gx, gy, z, ps = self._last_raster
        try:
            _scrivi_geotiff(fname, gx, gy, z, ps, self.epsg)
            QMessageBox.information(self, "Salvato", f"GeoTIFF salvato:\n{fname}")
        except Exception as e:
            QMessageBox.warning(self, "Errore", str(e))


# ============================================================
# SALVATAGGIO GEOTIFF
# ============================================================
def _scrivi_geotiff(path, grid_x, grid_y, z, pixel_size, epsg):
    if not HAS_GDAL:
        raise ImportError("GDAL non disponibile")
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(path, len(grid_x), len(grid_y), 1, gdal.GDT_Float32)
    ds.SetGeoTransform((grid_x.min(), pixel_size, 0,
                        grid_y.max(), 0, -pixel_size))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    band.WriteArray(np.flipud(z))
    band.SetNoDataValue(float('nan'))
    band.FlushCache()
    ds = None


# ============================================================
# INTEGRAZIONE QGIS
# ============================================================
def avvia_da_qgis():
    """
    Punto di ingresso quando lo script è lanciato dalla console QGIS.
    Legge il layer attivo e i punti selezionati, poi chiede il campo
    da interpolare e apre la finestra interattiva.
    """
    from qgis.PyQt.QtCore import QVariant

    layer = iface.activeLayer()
    if not layer:
        QMessageBox.warning(None, 'InSAR Suite',
                            'Nessun layer PS attivo.')
        return

    feats = list(layer.selectedFeatures())
    if not feats:
        QMessageBox.warning(None, 'InSAR Suite',
                            'Nessun punto PS selezionato.')
        return

    if len(feats) < MIN_PUNTI:
        QMessageBox.warning(
            None, 'InSAR Suite – Punti insufficienti',
            f'Seleziona almeno {MIN_PUNTI} punti PS per il kriging.\n'
            f'Punti selezionati: {len(feats)}'
        )
        return

    # ── Scelta campo da interpolare ───────────────────────────────────
    dlg_campo = SceltaCampoDialog()
    dlg_campo.populate(layer)
    if not dlg_campo.exec_():
        return
    scelta_tipo, scelta_campo = dlg_campo.getChoice()

    campi_d = [f.name() for f in layer.fields()
               if re.match(r'^D\d{8}$', f.name())]

    if scelta_tipo == 'ts' and not campi_d:
        QMessageBox.warning(None, 'InSAR Suite',
            'Nessun campo DYYYYMMDD trovato nel layer.\n'
            'I campi delle date devono avere formato DYYYYMMDD (es. D20170101).')
        return
    if scelta_tipo == 'field' and not scelta_campo:
        QMessageBox.warning(None, 'InSAR Suite', 'Nessun campo numerico selezionato.')
        return

    # Coordinate
    coords = np.array([[f.geometry().asPoint().x(),
                        f.geometry().asPoint().y()] for f in feats])

    # ── Calcolo valori ────────────────────────────────────────────────
    if scelta_tipo == 'ts':
        df = pd.DataFrame(
            [[f.id()] + [f[c] for c in campi_d] for f in feats],
            columns=['ID'] + campi_d
        )
        vals = df[campi_d].apply(pd.to_numeric, errors='coerce')
        t = np.array([
            (pd.to_datetime(c[1:], format='%Y%m%d') -
             pd.to_datetime(campi_d[0][1:], format='%Y%m%d')).days / 365.25
            for c in campi_d
        ])
        vel = []
        for i in range(len(df)):
            y = vals.iloc[i].values
            m = ~np.isnan(y)
            if np.sum(m) >= 2:
                a, _ = np.polyfit(t[m], y[m], 1)
                vel.append(a)
            else:
                vel.append(np.nan)
        vel = np.array(vel)
        label_val = 'Velocità TS (mm/anno)'
    else:
        vel = np.array([
            float(f[scelta_campo])
            if f[scelta_campo] is not None else np.nan
            for f in feats
        ])
        label_val = scelta_campo

    # ── Conversione coordinate ────────────────────────────────────────
    crs = layer.crs()
    if crs.isGeographic() and HAS_PYPROJ:
        lon, lat = coords[:, 0], coords[:, 1]
        zone = int(np.floor((np.mean(lon) + 180) / 6)) + 1
        epsg_code = 32600 + zone if np.mean(lat) >= 0 else 32700 + zone
        tr = Transformer.from_crs(crs.authid(),
                                   f'EPSG:{epsg_code}', always_xy=True)
        x, y = tr.transform(lon, lat)
        coords_m = np.column_stack((x, y))
    else:
        coords_m = coords
        epsg_code = crs.postgisSrid()

    dlg = KrigingInterattivoDialog(
        coords_m, vel, epsg=epsg_code,
        label=label_val,
        parent=iface.mainWindow()
    )
    dlg.exec_()


# ============================================================
# AVVIO
# ============================================================
avvia_da_qgis()
