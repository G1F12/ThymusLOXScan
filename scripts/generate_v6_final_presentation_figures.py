"""Create the three presentation-only v6 figures from frozen tabular outputs.

This script performs no expression-matrix access, statistical testing, barcode
matching, or scientific-table regeneration.  The only derived display value is
log2(CPM + 1) in Figure 2, calculated from the frozen per-sample CPM column.
"""

from __future__ import annotations

import csv
import io
import math
import subprocess
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "figures" / "v6_final"
FONT = Path(r"C:\Windows\Fonts\arial.ttf")
FONT_BOLD = Path(r"C:\Windows\Fonts\arialbd.ttf")
PHASE2_SOURCE = (
    "11ff0c7205bebb5f523bc0163bde6c2f59cef0d7:"
    "results/tables/gse240016_focal_raw_sample_summaries_v6.tsv"
)

NAVY = "#17365D"
BLUE = "#3C78A8"
PALE_BLUE = "#DCEAF5"
ORANGE = "#D07A2D"
PALE_ORANGE = "#F6E2CF"
GREEN = "#4F8A5B"
PALE_GREEN = "#DDEBDD"
GRAY = "#666666"
LIGHT = "#E8E8E8"
VERY_LIGHT = "#F7F7F7"
BLACK = "#171717"
WHITE = "#FFFFFF"


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(FONT_BOLD if bold else FONT), size)


def text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[float, float],
    value: str,
    size: int,
    *,
    bold: bool = False,
    fill: str = BLACK,
    anchor: str | None = None,
) -> None:
    draw.text(xy, value, font=font(size, bold), fill=fill, anchor=anchor)


def save(image: Image.Image, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    image.save(OUT / f"{stem}.png", dpi=(300, 300), optimize=True)
    image.convert("RGB").save(
        OUT / f"{stem}.pdf", "PDF", resolution=300.0, quality=95
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_phase2_sample_table() -> list[dict[str, str]]:
    raw = subprocess.check_output(
        ["git", "show", PHASE2_SOURCE], cwd=ROOT, text=True, encoding="utf-8"
    )
    return list(csv.DictReader(io.StringIO(raw), delimiter="\t"))


def figure_1() -> None:
    rows = read_tsv(ROOT / "release/v6.0-final/tables/final_focal_results_v6.tsv")
    expected = [
        ("capsFB", "Lox", -1.456, "3+3", "A"),
        ("medFB", "Loxl2", -0.992, "3+3", "B+"),
        ("medFB", "Loxl1", 0.752, "3+3", "B"),
        ("mTEC1", "Loxl2", -3.286, "2+2", "C"),
    ]
    observed = []
    for row in rows:
        subtype, gene = row["focal_comparison"].split()
        observed.append(
            (
                subtype,
                gene,
                round(float(row["primary_log2fc"]), 3),
                row["biological_n"],
                row["internal_grade"],
            )
        )
    if observed != expected:
        raise RuntimeError(f"Frozen focal table mismatch: {observed!r}")

    im = Image.new("RGB", (2400, 1450), WHITE)
    d = ImageDraw.Draw(im)
    text(d, (120, 85), "Four focal GSE240016 transcript estimates", 54, bold=True)
    text(
        d,
        (120, 155),
        "Aged versus young sample-level pseudobulk estimates; no confidence intervals were available.",
        30,
        fill=GRAY,
    )

    left, right, top, bottom = 650, 2220, 300, 1120
    xmin, xmax = -3.6, 1.1

    def xmap(value: float) -> float:
        return left + (value - xmin) / (xmax - xmin) * (right - left)

    for tick in [-3, -2, -1, 0, 1]:
        x = xmap(tick)
        d.line((x, top, x, bottom), fill=NAVY if tick == 0 else LIGHT, width=5 if tick == 0 else 2)
        text(d, (x, bottom + 40), f"{tick:+d}", 29, fill=GRAY, anchor="ma")
    text(d, ((left + right) / 2, bottom + 100), "Aged versus young log2 fold change", 34, bold=True, anchor="ma")

    y_positions = [410, 600, 790, 980]
    for (subtype, gene, value, n_value, grade), y in zip(observed, y_positions):
        label = f"{subtype}  {gene}"
        text(d, (120, y - 12), label, 38, bold=True, anchor="lm")
        text(d, (120, y + 42), f"n={n_value} sample-label proxies", 27, fill=GRAY, anchor="lm")
        x = xmap(value)
        color = BLUE if value < 0 else ORANGE
        d.line((xmap(0), y, x, y), fill=color, width=9)
        d.ellipse((x - 18, y - 18, x + 18, y + 18), fill=color, outline=WHITE, width=4)
        anchor = "rm" if value < 0 else "lm"
        offset = -28 if value < 0 else 28
        text(d, (x + offset, y - 42), f"{value:+.3f}", 32, bold=True, fill=color, anchor=anchor)
        text(d, (x + offset, y + 38), f"Category {grade}", 27, fill=GRAY, anchor=anchor)

    d.rounded_rectangle((120, 1260, 2280, 1370), radius=18, fill=VERY_LIGHT, outline=LIGHT, width=2)
    text(
        d,
        (160, 1315),
        "A, B+, B and C are internal computational robustness categories—not universal evidence grades.",
        29,
        fill=BLACK,
        anchor="lm",
    )
    save(im, "figure_1_focal_effects")


def figure_2() -> None:
    rows = [
        r
        for r in read_phase2_sample_table()
        if r["focal_comparison"] == "mTEC1 Loxl2"
        and r["included_as_biological_replicate"] == "YES"
    ]
    if len(rows) != 4 or [r["age"] for r in rows] != ["02mo", "02mo", "18mo", "18mo"]:
        raise RuntimeError("Expected exactly two young and two aged included sample proxies")

    im = Image.new("RGB", (2600, 1650), WHITE)
    d = ImageDraw.Draw(im)
    text(d, (120, 80), "mTEC1 Loxl2 sample-level architecture", 54, bold=True)
    text(
        d,
        (120, 150),
        "Four included sample-label proxies (2 young + 2 aged); cells are not independent replicates.",
        30,
        fill=GRAY,
    )

    labels = [r["sample"].replace("_d0", "") for r in rows]
    colors = [BLUE, BLUE, ORANGE, ORANGE]
    y_positions = [450, 690, 990, 1230]
    for y, r, label, color in zip(y_positions, rows, labels, colors):
        age = "Young (2 mo)" if r["age"] == "02mo" else "Aged (18 mo)"
        text(d, (120, y - 28), label, 40, bold=True, anchor="lm")
        text(d, (120, y + 28), age, 34, fill=color, anchor="lm")
        text(d, (120, y + 74), "biological-unit proxy", 30, fill=GRAY, anchor="lm")

    panel_lefts = [770, 1680]
    panel_width = 760
    panel_titles = ["Detection fraction", "log2(CPM + 1)"]
    maxima = [0.10, 3.6]
    values = [
        [float(r["detection_fraction"]) for r in rows],
        [math.log2(float(r["cpm"]) + 1.0) for r in rows],
    ]

    for left, title, max_value, panel_values in zip(panel_lefts, panel_titles, maxima, values):
        right = left + panel_width
        text(d, ((left + right) / 2, 300), title, 36, bold=True, anchor="ma")
        for frac in [0, 0.25, 0.5, 0.75, 1.0]:
            x = left + frac * panel_width
            d.line((x, 370, x, 1370), fill=LIGHT, width=2)
            tick = max_value * frac
            label = f"{tick:.3f}" if max_value < 1 else f"{tick:.1f}"
            text(d, (x, 1410), label, 31, fill=GRAY, anchor="ma")
        for y, value, color in zip(y_positions, panel_values, colors):
            x = left + value / max_value * panel_width
            d.line((left, y, x, y), fill=PALE_BLUE if color == BLUE else PALE_ORANGE, width=18)
            d.ellipse((x - 18, y - 18, x + 18, y + 18), fill=color, outline=WHITE, width=4)
            label = f"{value:.4f}" if max_value < 1 else f"{value:.2f}"
            text(d, (x + 25, y), label, 34, bold=True, fill=color, anchor="lm")

    d.line((90, 845, 2510, 845), fill=NAVY, width=3)
    text(d, (2480, 815), "2 + 2 design", 34, bold=True, fill=NAVY, anchor="ra")
    text(
        d,
        (120, 1560),
        "Display transform only: log2(CPM+1) was calculated from the frozen per-sample CPM column.",
        31,
        fill=GRAY,
    )
    save(im, "figure_2_mtec1_loxl2_sample_context")


def figure_3() -> None:
    rows = read_tsv(ROOT / "release/v6.0-final/tables/final_external_evidence_v6.tsv")
    by_claim: dict[str, list[dict[str, str]]] = {}
    claim_ids = {
        "capsFB Lox aged-lower": "C002",
        "medFB Loxl1 aged-higher": "C003",
        "medFB Loxl2 aged-lower": "C004",
        "mTEC1 Loxl2 aged-lower": "C005",
    }
    for row in rows:
        by_claim.setdefault(claim_ids[row["claim"]], []).append(row)

    row_specs = [
        ("C002", "capsFB Lox", "Lower", "Supportive broad context", "Unavailable / not eligible"),
        ("C003", "medFB Loxl1", "Higher", "Opposite broad context", "Unavailable / not eligible"),
        ("C004", "medFB Loxl2", "Lower", "Supportive broad context", "Unavailable / not eligible"),
        ("C005", "mTEC1 Loxl2", "Lower", "Supportive broad context", "Mixed related-subtype context"),
    ]
    for claim, *_ in row_specs:
        if claim not in by_claim:
            raise RuntimeError(f"Missing external evidence for {claim}")

    im = Image.new("RGB", (3000, 1800), WHITE)
    d = ImageDraw.Draw(im)
    text(d, (110, 70), "External evidence context for the four focal claims", 54, bold=True)
    text(
        d,
        (110, 140),
        "Independent datasets provide broad or related-subtype context; none establishes exact-subtype replication.",
        30,
        fill=GRAY,
    )

    x_edges = [100, 700, 1180, 1900, 2600, 2910]
    headers = [
        "Focal claim",
        "Primary\nGSE240016",
        "GSE223049\nbroad context",
        "E-MTAB-8560\nrelated TEC context",
        "Exact subtype",
    ]
    for i, header in enumerate(headers):
        x0, x1 = x_edges[i], x_edges[i + 1]
        d.rectangle((x0, 280, x1, 430), fill=NAVY, outline=WHITE, width=3)
        lines = header.split("\n")
        for j, line in enumerate(lines):
            header_size = 34 if i == 4 else 40
            text(d, ((x0 + x1) / 2, 340 + j * 48 - (len(lines) - 1) * 24), line, header_size, bold=True, fill=WHITE, anchor="mm")

    y_edges = [430, 690, 950, 1210, 1470]
    for row_index, (claim, focal, primary, gse, emtab) in enumerate(row_specs):
        y0, y1 = y_edges[row_index], y_edges[row_index + 1]
        cells = [
            (f"{claim}\n{focal}", VERY_LIGHT, BLACK),
            (primary, PALE_BLUE if primary == "Lower" else PALE_ORANGE, BLUE if primary == "Lower" else ORANGE),
            (gse, PALE_GREEN if gse.startswith("Supportive") else PALE_ORANGE, GREEN if gse.startswith("Supportive") else ORANGE),
            (emtab, VERY_LIGHT if emtab.startswith("Unavailable") else "#EEE6F5", GRAY if emtab.startswith("Unavailable") else "#6F4A8E"),
            ("Not established", VERY_LIGHT, GRAY),
        ]
        for i, (value, fill, color) in enumerate(cells):
            x0, x1 = x_edges[i], x_edges[i + 1]
            d.rectangle((x0, y0, x1, y1), fill=fill, outline=WHITE, width=5)
            lines = value.split("\n")
            cell_size = 44 if i == 0 else 36 if i == 4 else 42
            for j, line in enumerate(lines):
                text(
                    d,
                    ((x0 + x1) / 2, (y0 + y1) / 2 + j * 42 - (len(lines) - 1) * 21),
                    line,
                    cell_size,
                    bold=(i in {0, 1}),
                    fill=color,
                    anchor="mm",
                )

    d.rounded_rectangle((100, 1570, 2910, 1710), radius=18, fill=VERY_LIGHT, outline=LIGHT, width=2)
    text(d, (145, 1610), "GSE231906", 38, bold=True, fill=NAVY)
    text(d, (470, 1610), "DETECTION_CONTEXT_ONLY", 38, bold=True, fill=BLACK)
    text(
        d,
        (145, 1665),
        "Shown as a scope annotation only; it is excluded from mouse age-effect evidence.",
        36,
        fill=GRAY,
    )
    save(im, "figure_3_external_evidence_context")


if __name__ == "__main__":
    figure_1()
    figure_2()
    figure_3()
    print(f"Created six figure files under {OUT}")
