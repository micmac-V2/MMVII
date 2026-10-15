"""Import a GCP text file (the book's Reference ▸ Import Reference dialog)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog,
                               QFormLayout, QHBoxLayout, QHeaderView, QLabel, QLineEdit,
                               QPushButton, QSpinBox, QTableWidget, QTableWidgetItem, QVBoxLayout)

from .. import gcp as gcpmod
from .crs_widget import CRSEdit

DELIMITERS = [("Auto", "auto"), ("Comma", ","), ("Semicolon", ";"), ("Tab", "\t"), ("Whitespace", None)]
ROLE_LABELS = ["—", "name", "x", "y", "z", "sigma"]
ROLE_TIPS = "x = Easting or Longitude, y = Northing or Latitude, z = height, sigma = accuracy"
PREVIEW_ROWS = 30


class GCPImportDialog(QDialog):
    def __init__(self, project_crs: str, existing_sets: list[str], parent=None, path: str = ""):
        super().__init__(parent)
        self.setWindowTitle("Import ground control points")
        self.resize(820, 600)
        self._text = ""
        self._existing = set(existing_sets)
        self.points: list[gcpmod.GCP] = []
        self.warnings: list[str] = []

        self.path_edit = QLineEdit(path)
        browse = QPushButton("Browse…")
        browse.clicked.connect(self._browse)
        row = QHBoxLayout()
        row.addWidget(self.path_edit, 1)
        row.addWidget(browse)

        self.delim = QComboBox()
        for label, _ in DELIMITERS:
            self.delim.addItem(label)
        self.skip = QSpinBox()
        self.skip.setRange(0, 1000)
        self.skip.setToolTip("Header lines to skip (the book's 'Start import at row 2' = skip 1)")
        self.crs_edit = CRSEdit(project_crs)
        self.set_name = QLineEdit(self._default_set_name())
        self.sigma = QDoubleSpinBox()
        self.sigma.setRange(0.0001, 100.0)
        self.sigma.setDecimals(4)
        self.sigma.setValue(0.02)
        self.sigma.setSuffix(" (units of CRS)")
        self.sigma.setToolTip("Accuracy of each GCP coordinate, used when the file has no sigma column")
        self.project_crs = project_crs

        form = QFormLayout()
        form.addRow("File", row)
        form.addRow("Delimiter", self.delim)
        form.addRow("Skip header rows", self.skip)
        form.addRow("GCP coordinate system", self.crs_edit)
        form.addRow("MMVII set name", self.set_name)
        form.addRow("Default accuracy", self.sigma)

        self.table = QTableWidget()
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        self.status = QLabel()
        self.status.setWordWrap(True)

        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Import with MMVII")
        self.buttons.accepted.connect(self._accept)
        self.buttons.rejected.connect(self.reject)

        lay = QVBoxLayout(self)
        lay.addLayout(form)
        lay.addWidget(QLabel("Assign a role to each column with the boxes in the first row "
                             "(x = Easting/Longitude, y = Northing/Latitude):"))
        lay.addWidget(self.table, 1)
        lay.addWidget(self.status)
        lay.addWidget(self.buttons)

        self.path_edit.editingFinished.connect(self._load)
        self.delim.currentIndexChanged.connect(lambda _: self._refresh())
        self.skip.valueChanged.connect(self._validate)
        self.crs_edit.changed.connect(self._validate)
        self.set_name.textChanged.connect(self._validate)
        self._role_boxes: list[QComboBox] = []
        if path:
            self._load()
        else:
            self._validate()

    def _default_set_name(self) -> str:
        name, k = "GCP", 1
        while name in self._existing:
            k += 1
            name = f"GCP{k}"
        return name

    def _browse(self):
        f, _ = QFileDialog.getOpenFileName(self, "GCP file", self.path_edit.text(),
                                           "Text files (*.csv *.txt *.CSV *.TXT);;All files (*)")
        if f:
            self.path_edit.setText(f)
            self._load()

    def _load(self):
        p = Path(self.path_edit.text())
        if not p.is_file():
            self._text = ""
            self._refresh()
            return
        self._text = p.read_text(errors="replace")
        spec = gcpmod.guess_spec(self._text)
        self._guessed = spec
        idx = next((i for i, (_, d) in enumerate(DELIMITERS) if d == spec.delimiter and i > 0), 0)
        self.delim.blockSignals(True)
        self.delim.setCurrentIndex(idx)
        self.delim.blockSignals(False)
        self.skip.setValue(spec.skip_rows)
        self._refresh(spec)

    def _delimiter(self):
        d = DELIMITERS[self.delim.currentIndex()][1]
        return gcpmod.sniff_delimiter(self._text) if d == "auto" else d

    def _refresh(self, spec: gcpmod.TableSpec | None = None):
        rows = gcpmod.read_rows(self._text, self._delimiter()) if self._text else []
        ncol = max((len(r) for r in rows[:PREVIEW_ROWS]), default=0)
        old_roles = {k: b.currentText() for k, b in enumerate(self._role_boxes)}
        self.table.clear()
        self.table.setColumnCount(ncol)
        self.table.setRowCount(min(len(rows), PREVIEW_ROWS) + 1)
        self._role_boxes = []
        roles_by_col = {c: r for r, c in spec.columns.items()} if spec else {}
        for c in range(ncol):
            box = QComboBox()
            box.addItems(ROLE_LABELS)
            box.setToolTip(ROLE_TIPS)
            box.setCurrentText(roles_by_col.get(c, old_roles.get(c, "—")))
            box.currentIndexChanged.connect(self._validate)
            self.table.setCellWidget(0, c, box)
            self._role_boxes.append(box)
        for r, row in enumerate(rows[:PREVIEW_ROWS], start=1):
            for c, val in enumerate(row):
                self.table.setItem(r, c, QTableWidgetItem(val))
        self.table.setVerticalHeaderLabels(["role"] + [str(i) for i in range(1, self.table.rowCount())])
        self._validate()

    def spec(self) -> gcpmod.TableSpec:
        cols = {}
        for c, box in enumerate(self._role_boxes):
            role = box.currentText()
            if role != "—":
                cols[role] = c
        return gcpmod.TableSpec(delimiter=self._delimiter(), skip_rows=self.skip.value(), columns=cols)

    def _validate(self):
        ok_btn = self.buttons.button(QDialogButtonBox.StandardButton.Ok)
        msgs, ok = [], True
        roles = [b.currentText() for b in self._role_boxes if b.currentText() != "—"]
        if len(roles) != len(set(roles)):
            msgs.append("Each role can be used once.")
            ok = False
        if not self._text:
            msgs.append("Choose a file.")
            ok = False
        if not self.crs_edit.is_valid():
            msgs.append("Give the coordinate system the GCPs were surveyed in.")
            ok = False
        name = self.set_name.text().strip()
        if not name or not all(ch.isalnum() or ch in "_-" for ch in name):
            msgs.append("Set name: letters, digits, '_' or '-'.")
            ok = False
        elif name in self._existing:
            msgs.append(f"Set '{name}' exists and will be replaced.")
        if ok:
            try:
                self.points, self.warnings = gcpmod.parse_points(self._text, self.spec())
                msgs.insert(0, f"{len(self.points)} points read.")
                msgs += self.warnings[:5]
                if len(self.warnings) > 5:
                    msgs.append(f"… {len(self.warnings) - 5} more warnings")
                ok = bool(self.points)
            except ValueError as e:
                msgs.append(str(e))
                ok = False
        if ok and self.project_crs and self.crs_edit.text() != self.project_crs:
            msgs.append(f"MMVII will convert them to the project system {self.project_crs}.")
        self.status.setText("\n".join(msgs))
        ok_btn.setEnabled(ok)

    def _accept(self):
        self._validate()
        if self.buttons.button(QDialogButtonBox.StandardButton.Ok).isEnabled():
            self.accept()
