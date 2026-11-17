"""A line edit for a CRS (EPSG code, PROJ string, WKT) that shows the resolved name."""

from __future__ import annotations

from pyproj import CRS
from pyproj.exceptions import CRSError
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QLabel, QLineEdit, QVBoxLayout, QWidget

# Shortcuts shown as a tooltip; EPSG:6653 is the system used in the Saczuk manual.
EXAMPLES = ("EPSG:6653  NAD83(CSRS) / UTM 10N + CGVD2013 height\n"
            "EPSG:32610 WGS 84 / UTM 10N\n"
            "EPSG:4326  WGS 84 lat/lon\n"
            "IGNF:LAMB93, '+proj=utm +zone=31 +datum=WGS84'")


class CRSEdit(QWidget):
    changed = Signal()

    def __init__(self, value: str = "", placeholder: str = "e.g. EPSG:6653", parent=None):
        super().__init__(parent)
        self.edit = QLineEdit(value)
        self.edit.setPlaceholderText(placeholder)
        self.edit.setToolTip(EXAMPLES)
        self.label = QLabel()
        self.label.setWordWrap(True)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.edit)
        lay.addWidget(self.label)
        self.edit.textChanged.connect(self._update)
        self._update()

    def text(self) -> str:
        return self.edit.text().strip()

    def crs(self) -> CRS | None:
        try:
            return CRS.from_user_input(self.text()) if self.text() else None
        except CRSError:
            return None

    def is_valid(self) -> bool:
        return self.crs() is not None

    def _update(self):
        if not self.text():
            self.label.setText("<i>empty</i>")
        else:
            c = self.crs()
            self.label.setText(f"✔ {c.name}" if c else "<span style='color:#c00'>✘ not a known CRS</span>")
        self.changed.emit()
