import os
import re
import sys

from fontTools.ttLib import TTFont
from fontTools.pens.freetypePen import FreeTypePen
from fontTools.pens.boundsPen import BoundsPen
from fontTools.misc.transform import Transform

TARGET_WIDTH = 32
TARGET_HEIGHT = 64
REMOVE_DOT_RE = re.compile(r"\..*$")

def fit_to_higher_bitsize(n: int) -> int:
    if n <= 8:
        return 8
    if n <= 16:
        return 16
    if n <= 32:
        return 32
    if n <= 64:
        return 64
    raise ValueError("Value too big")


HEADER = f"""#include <array>
#include <map>
#include <cstdint>

namespace PlyFontData::[NAME] {"{"}

using LineType = uint{fit_to_higher_bitsize(TARGET_WIDTH)}_t;
static const int WIDTH = {TARGET_WIDTH};
static const int HEIGHT = {TARGET_HEIGHT};

static const std::map<unsigned char, std::array<LineType, HEIGHT>> RASTERIZED_CHARACTERS{"{"}"""

FOOTER = "};\n}\n"


def get_max_bounds(font):
    glyf_table = font["glyf"]
    glyph_set = font.getGlyphSet()
    cmap = font.getBestCmap()

    max_bounds = [0, 0, 0, 0]

    for glyph_name in cmap.values():
        if not glyf_table[glyph_name].isComposite():
            continue
        pen = BoundsPen(glyph_set)
        glyph_set[glyph_name].draw(pen)

        xmin, ymin, xmax, ymax = pen.bounds

        max_bounds[0] = min(max_bounds[0], xmin)
        max_bounds[1] = min(max_bounds[1], ymin)
        max_bounds[2] = max(max_bounds[2], xmax)
        max_bounds[3] = max(max_bounds[3], ymax)

    return tuple(max_bounds)


def code_gen(m: list[list[bool]], c: int):
    lines: list[int] = []

    for line in m:
        value = 0
        for i in range(len(line)):
            v = line[i]
            value <<= 1
            if not v:
                continue
            value |= 1
        lines.append(value)

    lines_hex: list[str] = [hex(value) for value in lines]

    lines_cpp = "{" + ", ".join(lines_hex) + "}"
    char_cpp = f"/* {chr(c)} */ {c}, "

    return "{" + char_cpp + lines_cpp + "}"


def byte_array_to_bool(arr: bytes, size: tuple[int, int]) -> list[list[bool]]:
    result: list[list[bool]] = []

    for y in range(size[1]):
        line: list[bool] = []
        for x in range(size[0]):
            line.append(arr[x + y * size[0]] > 127)
        result.append(line)

    return result


def get_matrix_by_ord(font_path: str) -> dict[int, list[list[bool]]]:
    result = {}

    with TTFont(font_path) as font:
        cmap = font.getBestCmap()

        glyph_set = font.getGlyphSet()

        xmin, ymin, xmax, ymax = get_max_bounds(font)

        w = xmax - xmin
        h = ymax - ymin

        ratio = min(TARGET_WIDTH / w, TARGET_HEIGHT / h)

        transform = Transform().scale(ratio, ratio).translate(0, h / 4)

        for glyph_ord, glyph_name in cmap.items():
            if glyph_ord > 255:
                continue

            pen = FreeTypePen(glyph_set)
            glyph_set[glyph_name].draw(pen)

            colors, size = pen.buffer(
                width=TARGET_WIDTH, height=TARGET_HEIGHT, transform=transform
            )
            result[glyph_ord] = byte_array_to_bool(colors, size)
    return result


def raster_to_stdout(raster: dict[int, list[list[bool]]]):
    for char, m in raster.items():
        print(f"--- {char} ({chr(char)}) ---")

        for line in m:
            for b in line:
                if b:
                    print("█", end="")
                else:
                    print(" ", end="")
            print()


def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print('Usage: generate.py <font file> <output file>')
        return

    font = args[0]
    output = args[1]
    debug = os.getenv("SHOW_FONT_RASTER") != None

    rasterized = get_matrix_by_ord(font)
    if debug:
        raster_to_stdout(rasterized)
        return
    
    with open(output, 'w') as f:
        gen_code(rasterized, f)


def gen_code(raster: dict[int, list[list[bool]]], f=sys.stdout):
    entries = []

    name = "no_name" if f == sys.stdout else REMOVE_DOT_RE.sub("", os.path.basename(f.name))
    for char, m in raster.items():
        entries.append(code_gen(m, char))

    print(HEADER.replace("[NAME]", name), file=f)
    print(",\n".join(entries), file=f)
    print(FOOTER, file=f)


if __name__ == "__main__":
    main()
