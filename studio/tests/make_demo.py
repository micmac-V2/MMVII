"""Build a demo project: synthetic RedEdge-M flight + GCP CSV in NAD83(CSRS) UTM 10N.

    python tests/make_demo.py /path/to/demo
"""

import sys
from pathlib import Path

from pyproj import Transformer

sys.path.insert(0, str(Path(__file__).parent))
import synth  # noqa: E402


def main(root: Path, n_captures: int = 36):
    synth.make_micasense(root, n_captures=n_captures)
    tr = Transformer.from_crs("EPSG:4326", "EPSG:3157", always_xy=True)
    pos = synth.flight_positions(n_captures)
    lines = ["Label,Easting,Northing,Height"]
    for k, i in enumerate((0, 7, 14, 21, 28, 35)[: max(1, n_captures // 6)], start=1):
        lat, lon = pos[min(i, n_captures - 1)]
        x, y = tr.transform(lon + 0.00007, lat - 0.00005)
        lines.append(f"{k},{x:.3f},{y:.3f},{4.5 + 0.3 * k:.3f}")
    (root / "Boundary.csv").write_text("\n".join(lines) + "\n")
    return root


if __name__ == "__main__":
    main(Path(sys.argv[1]))
