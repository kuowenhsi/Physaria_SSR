from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import rasterio
from rasterio.windows import from_bounds
from rasterio.warp import transform_bounds


INPUT_TABLE = Path(r"C:\Users\kuowe\Physaria_SSR\data\Table_1.csv")
INPUT_RASTER = Path(r"C:\Users\kuowe\Physaria_SSR\data\hfp_2020_100m_usa_extent.tif")
OUTPUT_RASTER = Path(r"C:\Users\kuowe\Physaria_SSR\data\hfp_2020_100m_p_globosa_water_overlay_extent.tif")


def parse_number(value: str) -> float:
    match = re.search(r"-?\d+(?:\.\d+)?", value)
    if not match:
        raise ValueError(f"Could not parse numeric value from {value!r}")
    return float(match.group(0))


def canonical_species(value: str) -> str:
    return " ".join(value.replace("\u201c", "").replace("\u201d", "").split())


def lonlat_to_web_mercator(lon: float, lat: float) -> tuple[float, float]:
    x = lon * 20037508.34 / 180.0
    y = math.log(math.tan((90.0 + lat) * math.pi / 360.0)) / (math.pi / 180.0)
    y = y * 20037508.34 / 180.0
    return x, y


def compute_overlay_bbox_3857() -> tuple[float, float, float, float]:
    xs: list[float] = []
    ys: list[float] = []

    with INPUT_TABLE.open(newline="", encoding="utf-8-sig") as handle:
        filtered_lines = []
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#") or stripped.startswith('"#'):
                continue
            filtered_lines.append(line)
        reader = csv.DictReader(filtered_lines)
        for row in reader:
            species = canonical_species(row["Species name"])
            if species != "P. globosa":
                continue
            lat = parse_number(row["Lat"])
            lon = parse_number(row["Long"])
            x, y = lonlat_to_web_mercator(lon, lat)
            xs.append(x)
            ys.append(y)

    if not xs or not ys:
        raise RuntimeError("No P. globosa coordinates were found in Table_1.csv")

    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    x_span = xmax - xmin
    y_span = ymax - ymin
    x_pad = x_span * 0.24
    y_pad = y_span * 0.26
    return xmin - x_pad, ymin - y_pad, xmax + x_pad, ymax + y_pad


def main() -> None:
    bbox_3857 = compute_overlay_bbox_3857()

    with rasterio.open(INPUT_RASTER) as src:
        bbox_hfp = transform_bounds("EPSG:3857", src.crs, *bbox_3857, densify_pts=100)
        window = from_bounds(*bbox_hfp, transform=src.transform)
        window = window.round_offsets().round_lengths()
        data = src.read(window=window)
        transform = src.window_transform(window)
        meta = src.meta.copy()

    OUTPUT_RASTER.parent.mkdir(parents=True, exist_ok=True)
    meta.update(
        {
            "height": data.shape[1],
            "width": data.shape[2],
            "transform": transform,
            "compress": "LZW",
            "tiled": True,
            "bigtiff": "IF_SAFER",
        }
    )

    with rasterio.open(OUTPUT_RASTER, "w", **meta) as dst:
        dst.write(data)

    print("Overlay bbox EPSG:3857:", bbox_3857)
    print("Overlay bbox HFP CRS:", bbox_hfp)
    print(f"Wrote {OUTPUT_RASTER}")


if __name__ == "__main__":
    main()
