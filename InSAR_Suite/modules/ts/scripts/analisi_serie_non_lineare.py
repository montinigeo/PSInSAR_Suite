import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os
from qgis.core import QgsTask, QgsMessageLog, Qgis, QgsApplication
from PyQt5.QtWidgets import QFileDialog, QMessageBox, QInputDialog
import mplcursors
import pwlf
from scipy.stats import linregress
import re  # Import per regex

# Registro globale per prevenire garbage collection dei task attivi
_active_tasks = []

# ================= FUNZIONE DI CORRELAZIONE =================
def corr_valid(x, y):
    mask = ~np.isnan(x) & ~np.isnan(y)
    if np.sum(mask) < 5:
        return np.nan
    return np.corrcoef(x[mask], y[mask])[0, 1]


# ================= MAIN =================
def main():
    soglia_corr_default = 0.85

    layer = iface.activeLayer()
    if not layer:
        QMessageBox.warning(None, 'InSAR TS',
            'Nessun layer PS attivo.\n'
            'Seleziona un layer PS puntuale nel pannello Layer prima di avviare l\'analisi.')
        return
    selected_features = layer.selectedFeatures()
    if not selected_features:
        QMessageBox.warning(None, 'InSAR TS – Nessun PS selezionato!',
            'Nessun punto PS selezionato nel layer attivo.\n\n'
            'Seleziona uno o più punti PS sulla mappa con gli strumenti di selezione di QGIS, '
            'poi avvia nuovamente l\'analisi.')
        return
    num_selected = len(selected_features)

    if num_selected == 1:
        soglia_corr = soglia_corr_default
    else:
        soglia_input, ok = QInputDialog.getDouble(
            None,
            "Soglia di correlazione",
            "Inserisci la soglia di correlazione (0\u20131):",
            soglia_corr_default, 0.0, 1.0, 2
        )
        soglia_corr = soglia_input if ok else soglia_corr_default

    campi_date = [f.name() for f in layer.fields() if re.match(r"^D\d{8}$", f.name())]
    if not campi_date:
        QMessageBox.warning(None, 'InSAR TS',
            'Nessun campo data trovato nel layer.\n'
            'I campi delle date devono avere formato DYYYYMMDD (es. D20170101).')
        return
    date = [pd.to_datetime(c[1:], format="%Y%m%d") for c in campi_date]

    records = []
    for feat in selected_features:
        code = feat["CODE"] if "CODE" in feat.fields().names() else feat.id()
        values = [feat[c] for c in campi_date]
        records.append([code] + values)
    df = pd.DataFrame(records, columns=["CODE"] + campi_date)

    task = AnalisiCinematicaTask(
        "InSAR TS - Analisi serie non lineare",
        df, date, soglia_corr, campi_date,
        salva_excel=False, percorso_excel=None
    )
    _active_tasks.append(task)  # previene garbage collection
    QgsApplication.taskManager().addTask(task)


# ================= TASK QGIS =================
# ================= TASK QGIS =================
class AnalisiCinematicaTask(QgsTask):
    def __init__(self, description, df, date, soglia_corr, campi_date, salva_excel, percorso_excel):
        super().__init__(description, QgsTask.CanCancel)
        self.df = df.copy()
        self.date = date
        self.soglia_corr = soglia_corr
        self.campi_date = campi_date
        self.salva_excel = salva_excel
        self.percorso_excel = percorso_excel
        self.result = None

    def run(self):
        try:
            valori = self.df[self.campi_date].apply(pd.to_numeric, errors='coerce')
            n = len(self.df)

            if n == 1:
                ps_coerenti = self.df.copy()
                corr_df = None
                msg_info = "ℹ️ Analisi di un singolo PS."
            else:
                corr_matrix = np.full((n, n), np.nan)
                for i in range(n):
                    serie_i = valori.iloc[i].values.astype(float)
                    for j in range(i, n):
                        serie_j = valori.iloc[j].values.astype(float)
                        c = corr_valid(serie_i, serie_j)
                        corr_matrix[i, j] = corr_matrix[j, i] = c

                corr_df = pd.DataFrame(corr_matrix, columns=self.df["CODE"], index=self.df["CODE"])
                mask_valid = (corr_df >= self.soglia_corr)
                coerenti = mask_valid.sum(axis=1) >= (n / 2)
                ps_coerenti = self.df.loc[coerenti.values].reset_index(drop=True)

                if len(ps_coerenti) == 0:
                    msg_info = f"⚠️ Nessun PS coerente trovato tra {n} punti selezionati."
                    self.result = (None, None, None, None, None, None, msg_info, Qgis.Warning, False)
                    return True
                elif len(ps_coerenti) == 1:
                    msg_info = f"ℹ️ Solo 1 PS coerente trovato su {n} selezionati."
                else:
                    msg_info = f"✅ Analisi completata: trovati {len(ps_coerenti)} PS coerenti su {n} selezionati."

            serie_coerenti = ps_coerenti[self.campi_date].to_numpy(dtype=float)
            serie_media = np.nanmean(serie_coerenti, axis=0)
            df_media = pd.DataFrame({
                "data": self.date,
                "deformazione_media": serie_media
            }).dropna().reset_index(drop=True)

            # --- FIT PIECEWISE LINEARE AUTOMATICO ---
            x = mdates.date2num(df_media["data"].values)
            y = df_media["deformazione_media"].values

            pwlf_model = pwlf.PiecewiseLinFit(x, y)
            max_segments = 5
            res_bic = []
            for i in range(2, max_segments + 1):
                try:
                    pwlf_model.fit(i)
                    rss = pwlf_model.rss
                    n_points = len(x)
                    k = 2 * i
                    bic = n_points * np.log(rss / n_points) + k * np.log(n_points)
                    res_bic.append((i, bic))
                except:
                    res_bic.append((i, np.inf))

            best_segments = min(res_bic, key=lambda t: t[1])[0] if res_bic else 2
            pwlf_model.fit(best_segments)
            breaks = pwlf_model.fit_breaks
            slopes = pwlf_model.slopes
            intercepts = pwlf_model.intercepts

            segmenti = []
            for i in range(best_segments):
                start_date = mdates.num2date(breaks[i]).strftime("%Y-%m-%d")
                end_date = mdates.num2date(breaks[i + 1]).strftime("%Y-%m-%d")
                segmenti.append({
                    "segmento": i + 1,
                    "data_inizio": start_date,
                    "data_fine": end_date,
                    "pendenza": slopes[i],
                    "intercetta": intercepts[i]
                })
            df_segmenti = pd.DataFrame(segmenti)

            # Salvataggio in Excel disabilitato
            self.salva_excel = False
            self.percorso_excel = None

            self.result = (ps_coerenti, df_media, df_segmenti, breaks, pwlf_model,
                           self.percorso_excel, msg_info, Qgis.Info, True)
            return True

        except Exception as e:
            QgsMessageLog.logMessage(f"Errore task: {str(e)}", "Cinematica", Qgis.Critical)
            return False

    def finished(self, result):
        if not result or self.result is None:
            QgsMessageLog.logMessage("❌ Task fallito", "Cinematica", Qgis.Critical)
            QMessageBox.critical(None, 'InSAR TS – Errore',
                'Elaborazione non completata. Controlla il log di QGIS per i dettagli.')
            return

        ps_coerenti, df_media, df_segmenti, breaks, pwlf_model, percorso_excel, msg_info, msg_level, do_plot = self.result

        QgsMessageLog.logMessage(msg_info, "Cinematica", msg_level)

        if percorso_excel:
            QgsMessageLog.logMessage(f"📁 File Excel salvato in: {percorso_excel}", "Cinematica", Qgis.Info)

        if not do_plot or df_media is None:
            QMessageBox.warning(None, 'InSAR TS – Nessun PS coerente trovato!',
                'Nessun PS coerente trovato tra i punti selezionati.\n\n'
                'Prova ad abbassare la soglia di correlazione oppure a selezionare '
                'un\'area con PS cinematicamente più omogenei.')
            return
        # --- PLOT ---
        plt.close('all')
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_facecolor("#f7f7f7")

        ax.plot(df_media["data"], df_media["deformazione_media"],
                'o', markersize=3, color='darkolivegreen', label='Serie storica media')

        x_num = mdates.date2num(df_media["data"])
        slope_tot, intercept_tot, r_value, p_value, std_err = linregress(
            x_num, df_media["deformazione_media"]
        )

        x_pred = np.linspace(mdates.date2num(df_media["data"].min()),
                             mdates.date2num(df_media["data"].max()), 500)
        dates_pred = mdates.num2date(x_pred)

        y_tot = slope_tot * x_pred + intercept_tot

        # --- RETTA TOTALE (verde tratteggiato sottile) ---
        ax.plot(dates_pred, y_tot, 'darkolivegreen', linewidth=1,
                label='Regressione lineare complessiva')

        # --- Piecewise ---
        y_pred = pwlf_model.predict(x_pred)
        ax.plot(dates_pred, y_pred, 'firebrick', linewidth=1.0,
                label=f'Tratti lineari individuati ({len(breaks) - 1} segmenti)')

        # Breakpoints
        ylim = ax.get_ylim()
        y_top = ylim[1]
        offset = (ylim[1] - ylim[0]) * 0.05
        for bp in breaks:
            bp_date = mdates.num2date(bp)
            ax.axvline(bp_date, color='grey', linestyle='--', alpha=0.7)
            ax.text(bp_date, y_top + offset, bp_date.strftime("%Y-%m-%d"),
                    rotation=90, verticalalignment='bottom',
                    horizontalalignment='center', fontsize=8, color='grey')

        # --- Box testo a destra ---
        box_pos = ax.get_position()
        ax.set_position([box_pos.x0, box_pos.y0, box_pos.width * 0.75, box_pos.height])

        ax_text = fig.add_axes([box_pos.x0 + box_pos.width * 0.78,
                                box_pos.y0, box_pos.width * 0.22, box_pos.height])
        ax_text.axis('off')

        ytxt = 1.0
        line_h = 0.04
        block_space = 0.10

        data_inizio_tot = df_media["data"].min().strftime("%Y-%m-%d")
        data_fine_tot = df_media["data"].max().strftime("%Y-%m-%d")

        # --- Regressione totale (verde) ---
        ax_text.text(0, ytxt, "Regressione lineare complessiva:",
                     fontsize=10, color='darkolivegreen', fontfamily='monospace')
        ytxt -= line_h

        ax_text.text(0, ytxt, f"  periodo: {data_inizio_tot} - {data_fine_tot}",
                     fontsize=10, color='darkolivegreen', fontfamily='monospace')
        ytxt -= line_h

        vel_tot_annua = round(slope_tot * 365.25, 2)
        r2_tot_rounded = round(r_value ** 2, 3)  # <-- qui calcolo r^2

        ax_text.text(0, ytxt, f"  velocità = {vel_tot_annua:.2f} mm/anno",
                     fontsize=10, color='darkolivegreen', fontfamily='monospace')
        ytxt -= line_h

        ax_text.text(0, ytxt, f"  r2 = {r2_tot_rounded:.3f}",
                     fontsize=10, color='darkolivegreen', fontfamily='monospace')  # <-- etichetta aggiornata
        ytxt -= block_space

        # --- Segmenti (rosso) ---
        for i, row in df_segmenti.iterrows():
            start_num = mdates.date2num(pd.to_datetime(row['data_inizio']))
            end_num = mdates.date2num(pd.to_datetime(row['data_fine']))

            mask_segmento = (x_num >= start_num) & (x_num <= end_num)

            if np.sum(mask_segmento) > 1:
                slope_seg, intercept_seg, r_seg, p_seg, std_err_seg = linregress(
                    x_num[mask_segmento], df_media["deformazione_media"].values[mask_segmento]
                )
            else:
                slope_seg, intercept_seg, r_seg = np.nan, np.nan, np.nan

            vel_annua = round(row['pendenza'] * 365.25, 2)
            r2_seg_rounded = round(r_seg ** 2, 3) if not np.isnan(r_seg) else np.nan  # <-- r^2 per segmento

            ax_text.text(0, ytxt, f"Segmento {row['segmento']}:",
                         fontsize=10, color='firebrick', fontfamily='monospace')
            ytxt -= line_h

            ax_text.text(0, ytxt, f"  {row['data_inizio']} - {row['data_fine']}",
                         fontsize=10, color='firebrick', fontfamily='monospace')
            ytxt -= line_h

            ax_text.text(0, ytxt, f"  velocità = {vel_annua:.2f} mm/anno",
                         fontsize=10, color='firebrick', fontfamily='monospace')
            ytxt -= line_h

            ax_text.text(0, ytxt, f"  r2 = {r2_seg_rounded if not np.isnan(r2_seg_rounded) else 'N/A'}",
                         fontsize=10, color='firebrick', fontfamily='monospace')
            ytxt -= block_space

        ax.set_xlabel("Data")
        ax.set_ylabel("Deformazione (mm)")
        ax.set_title("Analisi di linearità della serie storica esaminata")
        ax.grid(True)
        ax.legend()

        plt.show()



main()
