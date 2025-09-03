"""Geometry and baseline analysis utilities."""

from shapely.geometry import LineString, Polygon
from shapely.ops import split
from pyproj import Geod

geod = Geod(ellps="WGS84")


def baseline_length_within_polygon(baseline_coords, polygon_coords):
    """Return length of baseline segments lying inside the polygon in meters.

    Parameters
    ----------
    baseline_coords : list[tuple[float, float]]
        List of (lon, lat) pairs describing the baseline as a polyline.
    polygon_coords : list[tuple[float, float]]
        List of (lon, lat) pairs describing a closed polygon region.

    Returns
    -------
    float
        Total length of the portion of baseline inside the polygon, in meters.
    """
    baseline = LineString(baseline_coords)
    polygon = Polygon(polygon_coords)
    intersection = baseline.intersection(polygon)
    if intersection.is_empty:
        return 0.0
    # intersection may be multilinestring
    if intersection.geom_type == "LineString":
        return geod.geometry_length(intersection)
    elif intersection.geom_type == "MultiLineString":
        return sum(geod.geometry_length(geom) for geom in intersection.geoms)
    else:
        return 0.0


def split_baseline_by_polygon(baseline_coords, polygon_coords):
    """Return segments of baseline split by polygon boundary.

    Segments are categorized as inside or outside the polygon.
    """
    baseline = LineString(baseline_coords)
    polygon = Polygon(polygon_coords)
    # Split baseline by polygon boundary lines
    boundary = polygon.boundary
    pieces = split(baseline, boundary)
    segments = []
    for piece in pieces:
        inside = piece.within(polygon)
        segments.append((list(piece.coords), inside))
    return segments
