"""3D view of camera centres and GCPs (pyvista). Orientations come in phase 3."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QLabel, QVBoxLayout, QWidget

from .map_view import BAND_COLORS
from ..bands import band_sort_key


class SceneView(QWidget):
    """Created lazily: the VTK render window only starts when the tab is first shown."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._lay = QVBoxLayout(self)
        self._lay.setContentsMargins(0, 0, 0, 0)
        self.plotter = None
        self._pending = None
        self.note = QLabel()
        self.note.setWordWrap(True)
        self._lay.addWidget(self.note)

    def ensure_plotter(self):
        if self.plotter is not None:
            return
        try:
            from pyvistaqt import QtInteractor
            self.plotter = QtInteractor(self)
            self._lay.addWidget(self.plotter.interactor, 1)
        except Exception as e:  # no OpenGL context, e.g. headless
            self.note.setText(f"3D view unavailable: {e}")
            return
        if self._pending:
            self.set_data(*self._pending)

    def set_data(self, cams, gcps, note: str = ""):
        """cams: list of (band, x, y, z); gcps: list of (name, x, y, z)."""
        self.note.setText(note)
        self._pending = (cams, gcps, note)
        if self.plotter is None:
            return
        p = self.plotter
        p.clear()
        pts_all = []
        origin = None
        allxyz = [c[1:] for c in cams] + [g[1:] for g in gcps]
        if allxyz:
            origin = np.array(allxyz, dtype=float).mean(axis=0)
        for band in sorted({c[0] for c in cams}, key=band_sort_key):
            xyz = np.array([c[1:] for c in cams if c[0] == band], dtype=float) - origin
            pts_all.append(xyz)
            p.add_points(xyz, color=BAND_COLORS.get(band, "#555555"), point_size=8,
                         render_points_as_spheres=True, label=band)
        if gcps:
            xyz = np.array([g[1:] for g in gcps], dtype=float) - origin
            pts_all.append(xyz)
            p.add_points(xyz, color="#e8a400", point_size=14, render_points_as_spheres=True, label="GCP")
            p.add_point_labels(xyz, [g[0] for g in gcps], font_size=10, point_size=1,
                               shape_opacity=0.3, always_visible=True)
        if pts_all:
            p.add_legend(bcolor=None, face="circle")
            p.show_grid(xtitle="X - x0", ytitle="Y - y0", ztitle="Z - z0")
        p.reset_camera()
