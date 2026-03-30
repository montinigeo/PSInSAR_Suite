"""
runner_ewud.py  –  Decomposizione East-West / Up-Down da dati InSAR
Usa qgis:joinbylocationsummary (mean) invece di native:intersection +
statisticsbycategories: meno step, stesso risultato, molto più veloce.
Indici spaziali creati prima di ogni join per massima performance.
"""

import traceback, math
from qgis.PyQt.QtCore import QThread, pyqtSignal
from qgis.core import (
    QgsProcessingContext, QgsProcessingFeedback,
    QgsVectorLayer
)
import processing

PI = math.pi


class EwudRunner(QThread):
    progress = pyqtSignal(int)
    log      = pyqtSignal(str)
    finished = pyqtSignal(dict, list)
    error    = pyqtSignal(str)

    def __init__(self, params, parent=None):
        super().__init__(parent)
        self.params = params

    def _info(self, msg):
        self.log.emit(f'<span style="color:#aed6f1">&nbsp;&nbsp;{msg}</span>')

    def _step(self, msg, pct):
        self.log.emit(f'<span style="color:#5dade2">→ {msg}</span>')
        self.progress.emit(pct)

    def _run(self, alg, params, ctx, feedback):
        r = processing.run(alg, params,
                           context=ctx, feedback=feedback,
                           is_child_algorithm=False)
        return r['OUTPUT']

    def _index(self, layer, ctx, feedback):
        processing.run('native:createspatialindex', {'INPUT': layer},
                       context=ctx, feedback=feedback,
                       is_child_algorithm=False)

    def run(self):
        try:
            ctx      = QgsProcessingContext()
            feedback = _Feedback(self.progress, self.log, total_steps=18)
            p        = self.params

            griglia    = p['griglia_ricamp']
            id_griglia = p['id_griglia']
            vel_asc    = p['vel_asc']
            vel_desc   = p['vel_desc']
            theta_asc  = p['offnadir_asc']
            theta_desc = p['offnadir_desc']
            alpha_asc  = p['azimut_asc']
            alpha_desc = p['azimut_desc']
            sgn        = -1.0 if p.get('right_looking', True) else 1.0

            # ── Coefficienti geometrici ───────────────────────────────────────
            feedback.next_step('Pre-calcolo coefficienti geometrici…')
            Ca_e = sgn * math.sin(theta_asc  * PI/180) * math.cos(alpha_asc  * PI/180)
            Ca_u =       math.cos(theta_asc  * PI/180)
            Cd_e = sgn * math.sin(theta_desc * PI/180) * math.cos(alpha_desc * PI/180)
            Cd_u =       math.cos(theta_desc * PI/180)
            A        = Cd_e**2/Cd_u + Cd_u
            B        = Ca_e**2/Ca_u + Ca_u
            D        = Cd_e/Cd_u - Ca_e/Ca_u
            Ce_ratio = Ca_e / Ca_u
            feedback.pushInfo(f'Ca_e={Ca_e:.8f}  Ca_u={Ca_u:.8f}')
            feedback.pushInfo(f'Cd_e={Cd_e:.8f}  Cd_u={Cd_u:.8f}')
            feedback.pushInfo(f'A={A:.8f}  B={B:.8f}  D={D:.8f}')

            # ── Riproiezione PS nel CRS della griglia (se necessario) ─────────
            from qgis.core import QgsCoordinateReferenceSystem, QgsProject
            grid_crs = griglia.crs()
            ps_asc_w  = p['ps_asc']
            ps_desc_w = p['ps_desc']
            if ps_asc_w.crs() != grid_crs:
                feedback.pushInfo(
                    f'Riproiezione PS ascending da {ps_asc_w.crs().authid()} '
                    f'a {grid_crs.authid()} per join spaziale…')
                r = processing.run('native:reprojectlayer', {
                    'INPUT': ps_asc_w, 'TARGET_CRS': grid_crs,
                    'OUTPUT': 'TEMPORARY_OUTPUT',
                }, context=ctx, feedback=feedback, is_child_algorithm=False)
                ps_asc_w = r['OUTPUT']
            if ps_desc_w.crs() != grid_crs:
                feedback.pushInfo(
                    f'Riproiezione PS descending da {ps_desc_w.crs().authid()} '
                    f'a {grid_crs.authid()} per join spaziale…')
                r = processing.run('native:reprojectlayer', {
                    'INPUT': ps_desc_w, 'TARGET_CRS': grid_crs,
                    'OUTPUT': 'TEMPORARY_OUTPUT',
                }, context=ctx, feedback=feedback, is_child_algorithm=False)
                ps_desc_w = r['OUTPUT']

            # ── Indici spaziali ───────────────────────────────────────────────
            feedback.next_step('Indici spaziali su PS e griglia…')
            self._index(ps_asc_w,  ctx, feedback)
            self._index(ps_desc_w, ctx, feedback)
            self._index(griglia,   ctx, feedback)

            # ══════════════════════════════════════════════════════════════════
            # FASE 1 — Media Va e Vd per cella con joinbylocationsummary
            #   Questo algoritmo è molto più veloce di intersection +
            #   statisticsbycategories perché usa direttamente l'indice spaziale
            #   e calcola la media in un'unica passata.
            # ══════════════════════════════════════════════════════════════════

            feedback.next_step('Media Va e conteggio Na per cella (join spaziale ascending)…')
            stat_asc = self._run('qgis:joinbylocationsummary', {
                'INPUT':          griglia,
                'JOIN':           ps_asc_w,
                'JOIN_FIELDS':    [vel_asc],
                'PREDICATE':      [0],          # intersects
                'SUMMARIES':      [0, 6],       # count=0, mean=6
                'DISCARD_NOMATCH': False,
                'OUTPUT':         'TEMPORARY_OUTPUT',
            }, ctx, feedback)

            feedback.next_step('Media Vd e conteggio Nd per cella (join spaziale descending)…')
            stat_desc = self._run('qgis:joinbylocationsummary', {
                'INPUT':          griglia,
                'JOIN':           ps_desc_w,
                'JOIN_FIELDS':    [vel_desc],
                'PREDICATE':      [0],
                'SUMMARIES':      [0, 6],       # count=0, mean=6
                'DISCARD_NOMATCH': False,
                'OUTPUT':         'TEMPORARY_OUTPUT',
            }, ctx, feedback)

            # I campi prodotti si chiamano vel_asc_count, vel_asc_mean, vel_desc_count, vel_desc_mean
            va_field   = vel_asc  + '_mean'
            na_field   = vel_asc  + '_count'
            vd_field   = vel_desc + '_mean'
            nd_field   = vel_desc + '_count'

            # ── Join Va + Vd + Na + Nd sulla griglia ─────────────────────────
            feedback.next_step('Join Va + Vd + Na + Nd per cella…')
            joined = self._run('native:joinattributestable', {
                'INPUT':              stat_asc,
                'INPUT_2':            stat_desc,
                'FIELD':              id_griglia,
                'FIELD_2':            id_griglia,
                'FIELDS_TO_COPY':     [vd_field, nd_field],
                'METHOD':             1,
                'DISCARD_NONMATCHING': True,
                'PREFIX':             'desc_',
                'OUTPUT':             'TEMPORARY_OUTPUT',
            }, ctx, feedback)

            # Rinomina in Va, Vd, Na, Nd
            feedback.next_step('Rinomina campi Va, Vd, Na, Nd…')
            tmp = self._run('native:renametablefield', {
                'INPUT': joined, 'FIELD': va_field,
                'NEW_NAME': 'Va', 'OUTPUT': 'TEMPORARY_OUTPUT',
            }, ctx, feedback)
            tmp = self._run('native:renametablefield', {
                'INPUT': tmp, 'FIELD': na_field,
                'NEW_NAME': 'Na', 'OUTPUT': 'TEMPORARY_OUTPUT',
            }, ctx, feedback)
            tmp = self._run('native:renametablefield', {
                'INPUT': tmp, 'FIELD': 'desc_' + vd_field,
                'NEW_NAME': 'Vd', 'OUTPUT': 'TEMPORARY_OUTPUT',
            }, ctx, feedback)
            tmp = self._run('native:renametablefield', {
                'INPUT': tmp, 'FIELD': 'desc_' + nd_field,
                'NEW_NAME': 'Nd', 'OUTPUT': 'TEMPORARY_OUTPUT',
            }, ctx, feedback)
            tmp = self._run('native:renametablefield', {
                'INPUT': tmp, 'FIELD': id_griglia,
                'NEW_NAME': 'ID_griglia', 'OUTPUT': 'TEMPORARY_OUTPUT',
            }, ctx, feedback)

            # Tieni solo i campi necessari
            feedback.next_step('Selezione campi ID, Va, Vd, Na, Nd…')
            prev = self._run('native:retainfields', {
                'INPUT':  tmp,
                'FIELDS': ['ID_griglia', 'Va', 'Vd', 'Na', 'Nd'],
                'OUTPUT': 'TEMPORARY_OUTPUT',
            }, ctx, feedback)

            # ══════════════════════════════════════════════════════════════════
            # FASE 2 — Calcolo E, U, vel, ang2, vprev
            # ══════════════════════════════════════════════════════════════════

            feedback.next_step('Calcolo E (East-West)…')
            prev = self._run('native:fieldcalculator', {
                'INPUT': prev, 'OUTPUT': 'TEMPORARY_OUTPUT',
                'FIELD_NAME': 'E', 'FIELD_TYPE': 0,
                'FIELD_LENGTH': 10, 'FIELD_PRECISION': 1,
                'FORMULA': (
                    f'CASE\n'
                    f'WHEN Va IS NOT NULL AND Vd IS NOT NULL\n'
                    f'THEN round(({A:.10f}*Vd - {B:.10f}*Va) / {D:.10f}, 1)\n'
                    f'WHEN Va IS NULL AND Vd IS NOT NULL\n'
                    f'THEN round({Cd_e:.10f}*Vd, 1)\n'
                    f'WHEN Va IS NOT NULL AND Vd IS NULL\n'
                    f'THEN round({Ca_e:.10f}*Va, 1)\n'
                    f'ELSE 9999\n'
                    f'END'
                ),
            }, ctx, feedback)

            feedback.next_step('Calcolo U (Up-Down)…')
            prev = self._run('native:fieldcalculator', {
                'INPUT': prev, 'OUTPUT': 'TEMPORARY_OUTPUT',
                'FIELD_NAME': 'U', 'FIELD_TYPE': 0,
                'FIELD_LENGTH': 10, 'FIELD_PRECISION': 1,
                'FORMULA': (
                    f'CASE\n'
                    f'WHEN Va IS NOT NULL AND Vd IS NOT NULL\n'
                    f'THEN round({-Ce_ratio:.10f}*E + {B:.10f}*Va, 1)\n'
                    f'WHEN Va IS NULL AND Vd IS NOT NULL\n'
                    f'THEN round({Cd_u:.10f}*Vd, 1)\n'
                    f'WHEN Va IS NOT NULL AND Vd IS NULL\n'
                    f'THEN round({Ca_u:.10f}*Va, 1)\n'
                    f'ELSE 9999\n'
                    f'END'
                ),
            }, ctx, feedback)

            feedback.next_step('Calcolo vel, ang2, vprev…')
            prev = self._run('native:fieldcalculator', {
                'INPUT': prev, 'OUTPUT': 'TEMPORARY_OUTPUT',
                'FIELD_NAME': 'vel', 'FIELD_TYPE': 0,
                'FIELD_LENGTH': 10, 'FIELD_PRECISION': 1,
                'FORMULA': 'round(sqrt(E*E + U*U), 1)',
            }, ctx, feedback)

            prev = self._run('native:fieldcalculator', {
                'INPUT': prev, 'OUTPUT': 'TEMPORARY_OUTPUT',
                'FIELD_NAME': 'ang2', 'FIELD_TYPE': 1,
                'FIELD_LENGTH': 20, 'FIELD_PRECISION': 10,
                'FORMULA': (
                    f'CASE\n'
                    f'WHEN E=0 AND U=0 THEN 9999\n'
                    f'WHEN atan2(E,U)*180/{PI} >= 0 THEN atan2(E,U)*180/{PI}\n'
                    f'ELSE atan2(E,U)*180/{PI} + 360\n'
                    f'END'
                ),
            }, ctx, feedback)

            out_param = p.get('Parametri_def', 'TEMPORARY_OUTPUT')
            if out_param == 'TEMPORARY_OUTPUT':
                out_param = 'memory:'
            prev = self._run('native:fieldcalculator', {
                'INPUT': prev, 'OUTPUT': out_param,
                'FIELD_NAME': 'vprev', 'FIELD_TYPE': 2,
                'FIELD_LENGTH': 0, 'FIELD_PRECISION': 0,
                'FORMULA': (
                    "CASE\n"
                    "WHEN E = 9999 THEN 'NO_DATA'\n"
                    "WHEN (Va >= -2 AND Va <= 2 AND Vd >= -2 AND Vd <= 2) OR vel <= 2 THEN 'STABLE'\n"
                    "WHEN (ang2 >= 0 AND ang2 <= 45) OR ang2 > 315 THEN 'UP'\n"
                    "WHEN  ang2 >  45 AND ang2 <= 135 THEN 'EST'\n"
                    "WHEN  ang2 > 135 AND ang2 <= 225 THEN 'DOWN'\n"
                    "WHEN  ang2 > 225 AND ang2 <= 315 THEN 'WEST'\n"
                    "ELSE 'NO_DATA'\n"
                    "END"
                ),
            }, ctx, feedback)
            param_def = prev

            # ── Poligoni_EWUD: griglia + Va, Vd, E, U, vel, ang2, vprev ─────
            feedback.next_step('Join EWUD su griglia (Poligoni_EWUD)…')
            out_poly = p.get('Poligoni_ewud', 'TEMPORARY_OUTPUT')
            if out_poly == 'TEMPORARY_OUTPUT':
                out_poly = 'memory:'
            poly = self._run('native:joinattributestable', {
                'INPUT':               griglia,
                'INPUT_2':             param_def,
                'FIELD':               id_griglia,
                'FIELD_2':             'ID_griglia',
                'FIELDS_TO_COPY':      ['Va', 'Vd', 'Na', 'Nd', 'E', 'U', 'vel', 'ang2', 'vprev'],
                'METHOD':              1,
                'DISCARD_NONMATCHING': True,
                'PREFIX':              '',
                'OUTPUT':              out_poly,
            }, ctx, feedback)

            # ── Centroidi_EWUD ────────────────────────────────────────────────
            feedback.next_step('Calcolo centroidi (Centroidi_EWUD)…')
            out_centr = p.get('Centroidi_ewud', 'TEMPORARY_OUTPUT')
            if out_centr == 'TEMPORARY_OUTPUT':
                out_centr = 'memory:'
            centr = self._run('native:centroids', {
                'ALL_PARTS': True,
                'INPUT':  poly,
                'OUTPUT': out_centr,
            }, ctx, feedback)

            def _to_layer(val, name):
                if isinstance(val, QgsVectorLayer):
                    val.setName(name); return val
                elif isinstance(val, str) and val:
                    lyr = QgsVectorLayer(val, name, 'ogr')
                    return lyr if lyr.isValid() else None
                return None

            results = {
                'Poligoni_ewud':  poly,
                'Centroidi_ewud': centr,
            }
            layers = [
                ('Centroidi_EWUD', _to_layer(centr, 'Centroidi_EWUD')),
                ('Poligoni_EWUD',  _to_layer(poly,  'Poligoni_EWUD')),
            ]
            self.finished.emit(results, layers)

        except Exception:
            self.error.emit(traceback.format_exc())


class _Feedback(QgsProcessingFeedback):
    def __init__(self, progress_signal, log_signal, total_steps):
        super().__init__()
        self._prog  = progress_signal
        self._log   = log_signal
        self._total = total_steps
        self._step  = 0
        self._base  = 0

    def next_step(self, label):
        self._step += 1
        self._base  = int((self._step - 1) / self._total * 100)
        self._log.emit(f'<span style="color:#5dade2">→ {label}</span>')
        self._prog.emit(self._base)

    def setProgress(self, p):
        self._prog.emit(min(self._base + int(p / self._total), 99))

    def pushInfo(self, info):
        self._log.emit(f'<span style="color:#aed6f1">&nbsp;&nbsp;{info}</span>')

    def pushWarning(self, w):
        self._log.emit(f'<span style="color:#f39c12">&nbsp;&nbsp;⚠ {w}</span>')

    def reportError(self, err, fatal=False):
        self._log.emit(f'<span style="color:#e74c3c">&nbsp;&nbsp;✖ {err}</span>')
