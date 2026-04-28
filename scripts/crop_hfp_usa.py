from __future__ import annotations

from contextlib import ExitStack
from pathlib import Path

import rasterio
from rasterio.merge import merge
from rasterio.warp import transform_bounds


SOURCE_DIR = Path(r"G:\GIS_Layers\HFP-100m-2020\2020")
OUTPUT_PATH = Path(r"C:\Users\kuowe\Physaria_SSR\data\hfp_2020_100m_usa_extent.tif")

# Bounding box for the 50 U.S. states in lon/lat, broad enough to include Alaska and Hawaii.
USA_BOUNDS_LONLAT = (-179.0, 18.8, -66.5, 71.5)


def intersects(a: tuple[float, float, float, float], b: tuple[float, float, float, float]) -> bool:
    return not (a[2] <= b[0] or a[0] >= b[2] or a[3] <= b[1] or a[1] >= b[3])


def main() -> None:
    tif_paths = sorted(SOURCE_DIR.glob("*.tif"))
    if not tif_paths:
        raise FileNotFoundError(f"No GeoTIFF tiles found in {SOURCE_DIR}")

    with rasterio.open(tif_paths[0]) as sample:
        target_bounds = transform_bounds(
            "EPSG:4326",
            sample.crs,
            *USA_BOUNDS_LONLAT,
            densify_pts=100,
        )
        selected_paths = []
        for path in tif_paths:
            with rasterio.open(path) as src:
                if intersects(tuple(src.bounds), target_bounds):
                    selected_paths.append(path)

        if not selected_paths:
            raise RuntimeError("No source tiles intersect the requested USA extent.")

        with ExitStack() as stack:
            datasets = [stack.enter_context(rasterio.open(path)) for path in selected_paths]
            mosaic, transform = merge(
                datasets,
                bounds=target_bounds,
                nodata=sample.nodata,
            )
            meta = sample.meta.copy()

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    meta.update(
        {
            "driver": "GTiff",
            "height": mosaic.shape[1],
            "width": mosaic.shape[2],
            "transform": transform,
            "compress": "LZW",
            "tiled": True,
            "bigtiff": "IF_SAFER",
        }
    )

    with rasterio.open(OUTPUT_PATH, "w", **meta) as dst:
        dst.write(mosaic)

    print(f"Selected {len(selected_paths)} tiles")
    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
