"""Main window: project folder, image table, GCP table, map and 3D tabs, MMVII log."""

from __future__ import annotations

from pathlib import Path

from pyproj import CRS, Transformer
from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (QApplication, QComboBox, QDialog, QDialogButtonBox, QDockWidget,
                               QFileDialog, QHBoxLayout, QHeaderView, QLabel, QMainWindow,
                               QMessageBox, QPlainTextEdit, QProgressDialog, QTableWidget,
                               QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from .. import gcp as gcpmod
from .. import mmvii
from ..project import Project, horizontal, missing_grids, reproject_xy, to_map_xy
from .crs_widget import CRSEdit
from .gcp_dialog import GCPImportDialog
from .map_view import MapCamera, MapGCP, MapView
from .scene_view import SceneView

ALL_BANDS = "All bands"


def _fmt(v, nd=3):
    return "" if v is None else f"{v:.{nd}f}"


def _table(headers: list[str]) -> QTableWidget:
    t = QTableWidget(0, len(headers))
    t.setHorizontalHeaderLabels(headers)
    t.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
    t.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
    t.horizontalHeader().setSortIndicator(0, Qt.SortOrder.AscendingOrder)
    t.setSortingEnabled(True)
    t.verticalHeader().setVisible(False)
    t.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
    return t


def _fill(table: QTableWidget, rows: list[list], keys: list[str] | None = None):
    table.setSortingEnabled(False)
    table.setRowCount(len(rows))
    for r, row in enumerate(rows):
        for c, val in enumerate(row):
            item = QTableWidgetItem()
            if isinstance(val, (int, float)) and not isinstance(val, bool):
                item.setData(Qt.ItemDataRole.DisplayRole, val)
            else:
                item.setText("" if val is None else str(val))
            if keys and c == 0:
                item.setData(Qt.ItemDataRole.UserRole, keys[r])
            table.setItem(r, c, item)
    table.setSortingEnabled(True)


class CRSDialog(QDialog):
    def __init__(self, current: str, auto_name: str, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Project coordinate system")
        self.edit = CRSEdit(current, placeholder=f"empty = automatic ({auto_name})")
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel("Coordinate system for the map, GCPs and later MMVII outputs.\n"
                             "Leave empty to use the UTM zone of the image GPS positions."))
        lay.addWidget(self.edit)
        lay.addWidget(bb)
        self.resize(480, 0)
        self.edit.changed.connect(lambda: bb.button(QDialogButtonBox.StandardButton.Ok)
                                  .setEnabled(not self.edit.text() or self.edit.is_valid()))


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("MMVII Studio")
        self.resize(1400, 900)
        self.project: Project | None = None
        self.settings = QSettings("MMVII", "Studio")

        # central tabs
        self.tabs = QTabWidget()
        self.map_view = MapView()
        self.map_view.imageClicked.connect(self.select_image)
        self.scene_view = SceneView()
        self.cameras_table = _table(["Make", "Model", "Band", "Width", "Height", "Focal (mm)",
                                     "35 mm equiv.", "Pixel (µm)", "Images"])
        self.tabs.addTab(self.map_view, "Map")
        self.tabs.addTab(self.scene_view, "3D")
        self.tabs.addTab(self.cameras_table, "Cameras")
        self.tabs.currentChanged.connect(self._tab_changed)
        self.setCentralWidget(self.tabs)

        # images dock (the book's Reference pane, cameras part)
        self.band_filter = QComboBox()
        self.band_filter.currentIndexChanged.connect(self.refresh_views)
        self.images_table = _table(["Image", "Band", "Lon", "Lat", "Alt (m)", "X", "Y", "Focal (mm)", "Time"])
        self.images_table.itemSelectionChanged.connect(self._image_selected)
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(4, 4, 4, 4)
        row = QHBoxLayout()
        row.addWidget(QLabel("Show"))
        row.addWidget(self.band_filter, 1)
        lay.addLayout(row)
        lay.addWidget(self.images_table)
        self.images_dock = self._dock("Images", w, Qt.DockWidgetArea.LeftDockWidgetArea)

        # GCP dock
        self.gcp_table = _table(["Set", "Name", "E / Lon", "N / Lat", "Z", "σ"])
        self.gcp_dock = self._dock("Ground control points", self.gcp_table, Qt.DockWidgetArea.LeftDockWidgetArea)

        # log dock
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(5000)
        self.log_dock = self._dock("MMVII log", self.log, Qt.DockWidgetArea.BottomDockWidgetArea)

        self._make_menus()
        self.status_label = QLabel()
        self.statusBar().addPermanentWidget(self.status_label)
        self._update_actions()
        exe = mmvii.find_mmvii()
        self.append_log(f"MMVII binary: {exe}" if exe else "MMVII binary not found: set MMVII_BIN.")

    # ---------------- UI construction ----------------
    def _dock(self, title, widget, area):
        d = QDockWidget(title, self)
        d.setObjectName(title)
        d.setWidget(widget)
        self.addDockWidget(area, d)
        return d

    def _make_menus(self):
        m = self.menuBar().addMenu("&Project")
        self.act_open = self._action(m, "&Open folder…", self.choose_folder, QKeySequence.StandardKey.Open)
        self.recent_menu = m.addMenu("Open &recent")
        self._fill_recent()
        m.addSeparator()
        self.act_crs = self._action(m, "Coordinate &system…", self.choose_crs)
        m.addSeparator()
        self._action(m, "&Quit", self.close, QKeySequence.StandardKey.Quit)

        m = self.menuBar().addMenu("&Import")
        self.act_scan = self._action(m, "&Scan images", self.scan_images, "Ctrl+I")
        self.act_gcp = self._action(m, "&Ground control points…", self.import_gcps, "Ctrl+G")

        m = self.menuBar().addMenu("&View")
        for d in (self.images_dock, self.gcp_dock, self.log_dock):
            m.addAction(d.toggleViewAction())

    def _action(self, menu, text, slot, shortcut=None):
        a = QAction(text, self)
        if shortcut:
            a.setShortcut(shortcut)
        a.triggered.connect(slot)
        menu.addAction(a)
        return a

    def _fill_recent(self):
        self.recent_menu.clear()
        for p in self.settings.value("recent", [], type=list):
            self.recent_menu.addAction(p, lambda p=p: self.open_project(Path(p)))
        self.recent_menu.setEnabled(not self.recent_menu.isEmpty())

    def _update_actions(self):
        has = self.project is not None
        for a in (self.act_crs, self.act_scan, self.act_gcp):
            a.setEnabled(has)

    def append_log(self, text: str):
        self.log.appendPlainText(text.rstrip())

    # ---------------- project ----------------
    def choose_folder(self):
        start = self.project.root if self.project else str(Path.home())
        d = QFileDialog.getExistingDirectory(self, "Project folder (contains the images)", str(start))
        if d:
            self.open_project(Path(d))

    def open_project(self, root: Path):
        if not root.is_dir():
            QMessageBox.warning(self, "Open project", f"{root} is not a folder.")
            return
        self.project = Project.open(root)
        recent = [str(root)] + [p for p in self.settings.value("recent", [], type=list) if p != str(root)]
        self.settings.setValue("recent", recent[:10])
        self._fill_recent()
        self.setWindowTitle(f"MMVII Studio — {root.name}")
        self.append_log(f"Opened project {root}")
        self._update_actions()
        if not self.project.images:
            self.scan_images()
        else:
            self.refresh_all()

    def choose_crs(self):
        p = self.project
        auto = Project(root=p.root, images=p.images).effective_crs().name
        dlg = CRSDialog(p.crs, auto, self)
        if dlg.exec():
            p.crs = dlg.edit.text()
            p.save()
            self.append_log(f"Project coordinate system: {p.crs or 'automatic'} ({p.effective_crs().name})")
            self.refresh_all()

    def scan_images(self):
        p = self.project
        dlg = QProgressDialog("Reading image metadata…", "Cancel", 0, 0, self)
        dlg.setWindowModality(Qt.WindowModality.WindowModal)
        dlg.setMinimumDuration(300)

        def progress(i, n, name):
            dlg.setMaximum(n)
            dlg.setValue(i)
            dlg.setLabelText(f"Reading {name} ({i + 1}/{n})")
            QApplication.processEvents()
            return not dlg.wasCanceled()

        p.scan_images(progress)
        dlg.close()
        p.save()
        n_gps = sum(im.has_gps for im in p.images)
        self.append_log(f"Scanned {len(p.images)} images in {len(p.captures)} captures, "
                        f"{n_gps} with GPS; bands: {', '.join(p.bands) or '-'}")
        for e in p.scan_errors[:20]:
            self.append_log(f"  skipped {e}")
        self.refresh_all()

    def import_gcps(self):
        p = self.project
        sets = [s.name for s in p.gcp_sets()]
        last = self.settings.value("lastGCPFile", "", type=str)
        dlg = GCPImportDialog(p.crs, sets, self, path=last)
        if not dlg.exec():
            return
        self.settings.setValue("lastGCPFile", dlg.path_edit.text())
        set_name = dlg.set_name.text().strip()
        try:
            QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
            out = gcpmod.import_with_mmvii(p.root, set_name, dlg.points, dlg.crs_edit.text(),
                                           p.crs or None, dlg.sigma.value())
            self.append_log(out)
        except mmvii.MMVIIError as e:
            self.append_log(e.output)
            QMessageBox.critical(self, "ImportGCP failed", f"{e}\n\nSee the MMVII log for details.")
        except (FileNotFoundError, ValueError) as e:
            QMessageBox.critical(self, "Import GCPs", str(e))
        finally:
            QApplication.restoreOverrideCursor()
        for w in dlg.warnings:
            self.append_log(f"  {w}")
        self.refresh_all()

    # ---------------- views ----------------
    def refresh_all(self):
        p = self.project
        cur = self.band_filter.currentText()
        self.band_filter.blockSignals(True)
        self.band_filter.clear()
        bands = p.bands if p else []
        self.band_filter.addItems(([ALL_BANDS] if len(bands) > 1 else []) + bands)
        # Default to one band: the book processes each band separately, starting with Blue.
        default = cur if cur in bands or cur == ALL_BANDS else (bands[0] if bands else "")
        self.band_filter.setCurrentText(default)
        self.band_filter.blockSignals(False)
        # "—" = not in the image metadata (e.g. early MicaSense firmware writes no Make/Model)
        _fill(self.cameras_table, [[v if v not in (None, "") else "—" for v in
                                    [c["make"], c["model"], c["band"], c["width"], c["height"],
                                     c["focal_mm"], c["focal_35mm"],
                                     round(c["pixel_size_mm"] * 1000, 2) if c["pixel_size_mm"] else None,
                                     c["count"]]] for c in (p.cameras_summary() if p else [])])
        self.refresh_views()

    def refresh_views(self):
        p = self.project
        if p is None:
            return
        crs = p.effective_crs()
        hcrs = horizontal(crs)
        units = hcrs.axis_info[0].unit_name if hcrs.axis_info else ""
        band = self.band_filter.currentText()
        imgs = [im for im in p.images if band in (ALL_BANDS, "") or im.band == band]
        gps = [im for im in imgs if im.has_gps]
        xy = to_map_xy([(im.lon, im.lat) for im in gps], crs) if gps else []
        pos = {im.path: (x, y) for im, (x, y) in zip(gps, xy)}

        _fill(self.images_table, [[Path(im.path).name, im.band, _fmt(im.lon, 7), _fmt(im.lat, 7),
                                   _fmt(im.alt), _fmt(pos.get(im.path, (None,))[0]),
                                   _fmt(pos.get(im.path, (None, None))[1]), im.focal_mm, im.datetime]
                                  for im in imgs], keys=[im.path for im in imgs])

        gcp_rows, map_gcps, scene_gcps = [], [], []
        for s in p.gcp_sets():
            scrs = crs
            try:
                scrs = CRS.from_user_input(s.crs) if s.crs else crs
                gxy = reproject_xy([(g.x, g.y) for g in s.points], scrs, crs)
            except Exception as e:
                self.append_log(f"GCP set {s.name}: cannot convert from '{s.crs}': {e}")
                gxy = [(g.x, g.y) for g in s.points]
            nd = 8 if scrs.is_geographic else 3
            for g, (x, y) in zip(s.points, gxy):
                gcp_rows.append([s.name, g.name, _fmt(g.x, nd), _fmt(g.y, nd), _fmt(g.z), _fmt(g.sigma, 4)])
                map_gcps.append(MapGCP(g.name, x, y))
                scene_gcps.append((g.name, x, y, g.z))
        _fill(self.gcp_table, gcp_rows)

        cams = [MapCamera(im.path, im.band, *pos[im.path]) for im in gps]
        self.map_view.set_data(cams, map_gcps, crs.name, f"({units})" if units else "")
        self._update_scene(gps, pos, crs, scene_gcps)
        n_gcp = len(map_gcps)
        self.status_label.setText(f"{p.root}  |  {crs.name}  |  {len(p.images)} images, "
                                  f"{len(p.captures)} captures  |  {n_gcp} GCPs")

    def _update_scene(self, gps, pos, crs, scene_gcps):
        note = ""
        heights = [im.alt or 0.0 for im in gps]
        if gps and crs.is_compound:
            src = CRS.from_epsg(4979)  # GPS: WGS 84 lat/lon + ellipsoidal height
            grids = missing_grids(src, crs)
            if grids:
                note = (f"Camera heights are raw GPS ellipsoidal heights, not heights in {crs.name}: PROJ is "
                        f"missing {', '.join(grids)}. Download it from https://cdn.proj.org into the "
                        "PROJ data folder (projsync --file NAME) — MMVII needs it too.")
            else:
                tr = Transformer.from_crs(src, crs, always_xy=True)
                heights = [tr.transform(im.lon, im.lat, h)[2] for im, h in zip(gps, heights)]
        cams = [(im.band, *pos[im.path], h) for im, h in zip(gps, heights)]
        self.scene_view.set_data(cams, scene_gcps, note)

    def _tab_changed(self, idx):
        if self.tabs.widget(idx) is self.scene_view:
            self.scene_view.ensure_plotter()

    def _image_selected(self):
        items = self.images_table.selectedItems()
        path = None
        if items:
            first = self.images_table.item(items[0].row(), 0)
            path = first.data(Qt.ItemDataRole.UserRole)
        self.map_view.select(path)

    def select_image(self, path: str):
        for r in range(self.images_table.rowCount()):
            if self.images_table.item(r, 0).data(Qt.ItemDataRole.UserRole) == path:
                self.images_table.selectRow(r)
                self.images_table.scrollToItem(self.images_table.item(r, 0))
                return
