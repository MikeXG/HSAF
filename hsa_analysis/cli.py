"""Command line interface for baseline analysis."""

import argparse
from .analysis import baseline_length_within_polygon


def parse_coord_pairs(values):
    if len(values) % 2 != 0:
        raise ValueError("Coordinate list must contain an even number of values")
    coords = []
    for i in range(0, len(values), 2):
        lon = float(values[i])
        lat = float(values[i + 1])
        coords.append((lon, lat))
    return coords


def main(argv=None):
    parser = argparse.ArgumentParser(description="Baseline and polygon intersection analysis")
    parser.add_argument("--baseline", nargs="+", type=float, help="Baseline lon lat pairs", required=True)
    parser.add_argument("--polygon", nargs="+", type=float, help="Polygon lon lat pairs", required=True)
    args = parser.parse_args(argv)

    baseline_coords = parse_coord_pairs(args.baseline)
    polygon_coords = parse_coord_pairs(args.polygon)
    length = baseline_length_within_polygon(baseline_coords, polygon_coords)
    print(f"Intersection length: {length:.2f} m")


if __name__ == "__main__":
    main()
