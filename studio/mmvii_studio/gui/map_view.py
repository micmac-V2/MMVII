"""2D map of camera centres and GCPs in the project CRS (the book's 'camera centres' view)."""

from __future__ import annotations

from dataclasses import dataclass

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QVBoxLayout, QWidget

from ..bands import band_sort_key

# Colours per band, matched to the band where it means something.
BAND_COLORS = {"RGB": "#4c72b0", "Blue": "#1f5fd6", "Green": "#2ca02c", "Red": "#d62728",
               "RedEdge": "#b5446e", "NIR": "#8c564b", "LWIR": "#ff7f0e", "Panchro": "#7f7f7f"}


@dataclass
class MapCamera:
    path: str
    band: str
    x: float
    y: float


@dataclass
class MapGCP:
    name: str
    x: float
    y: float


class MapView(QWidget):
    imageClicked = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(figsize=(6, 5), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.toolbar)
        lay.addWidget(self.canvas)
        self.ax = self.figure.add_subplot()
        self._cams: list[MapCamera] = []
        self._scatter_paths: dict[object, list[str]] = {}
        self._selected = None
        self.canvas.mpl_connect("pick_event", self._on_pick)
        self.set_data([], [], "", "")

    def set_data(self, cams: list[MapCamera], gcps: list[MapGCP], crs_name: str, units: str):
        self._cams = cams
        ax = self.ax
        ax.clear()
        self._scatter_paths.clear()
        self._selected = None
        bands = sorted({c.band for c in cams}, key=band_sort_key)
        for band in bands:
            sub = [c for c in cams if c.band == band]
            sc = ax.scatter([c.x for c in sub], [c.y for c in sub], s=18,
                            c=BAND_COLORS.get(band, "#555555"), label=f"{band} ({len(sub)})",
                            picker=5, zorder=2, edgecolors="none")
            self._scatter_paths[sc] = [c.path for c in sub]
        if len(cams) > 1:
            ax.plot([c.x for c in cams if c.band == bands[0]], [c.y for c in cams if c.band == bands[0]],
                    lw=0.5, c="#999999", zorder=1)
        if gcps:
            ax.scatter([g.x for g in gcps], [g.y for g in gcps], marker="^", s=70, c="#e8a400",
                       edgecolors="black", linewidths=0.6, label=f"GCP ({len(gcps)})", zorder=3)
            for g in gcps:
                ax.annotate(g.name, (g.x, g.y), xytext=(4, 4), textcoords="offset points", fontsize=8)
        if cams or gcps:
            ax.legend(loc="best", fontsize=8)
            ax.set_aspect("equal", adjustable="datalim")
        else:
            ax.text(0.5, 0.5, "No images with GPS yet.\nUse Project ▸ Open folder, then Import ▸ Scan images.",
                    ha="center", va="center", transform=ax.transAxes, color="#666666")
        ax.set_title(crs_name, fontsize=9)
        ax.set_xlabel(f"X / Easting {units}".strip())
        ax.set_ylabel(f"Y / Northing {units}".strip())
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.grid(True, lw=0.3)
        self.canvas.draw_idle()

    def select(self, path: str | None):
        if self._selected is not None:
            self._selected.remove()
            self._selected = None
        cam = next((c for c in self._cams if c.path == path), None)
        if cam:
            (self._selected,) = self.ax.plot([cam.x], [cam.y], marker="o", ms=12, mfc="none",
                                             mec="black", mew=1.5, zorder=4)
        self.canvas.draw_idle()

    def _on_pick(self, event):
        paths = self._scatter_paths.get(event.artist)
        if paths and len(event.ind):
            self.imageClicked.emit(paths[event.ind[0]])
