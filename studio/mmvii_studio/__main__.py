"""Entry point: `python -m mmvii_studio [project-folder]`."""

import sys
from pathlib import Path


def main():
    from PySide6.QtWidgets import QApplication

    from .gui.main_window import MainWindow

    app = QApplication(sys.argv)
    app.setApplicationName("MMVII Studio")
    win = MainWindow()
    win.show()
    if len(sys.argv) > 1:
        win.open_project(Path(sys.argv[1]))
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
