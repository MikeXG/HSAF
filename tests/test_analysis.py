import math
from shapely.geometry import LineString
from hsa_analysis.analysis import baseline_length_within_polygon, split_baseline_by_polygon


def test_baseline_length_within_polygon():
    baseline = [(110.0, 10.0), (112.0, 10.0), (112.0, 12.0)]
    polygon = [(111.0, 9.0), (113.0, 9.0), (113.0, 11.0), (111.0, 11.0)]
    length = baseline_length_within_polygon(baseline, polygon)
    # length should be approximately 222 km (2 degrees of longitude) along the equator
    assert math.isclose(length / 1000, 222, rel_tol=0.05)


def test_split_baseline_by_polygon():
    baseline = [(110.0, 10.0), (112.0, 10.0), (112.0, 12.0)]
    polygon = [(111.0, 9.0), (113.0, 9.0), (113.0, 11.0), (111.0, 11.0)]
    segments = split_baseline_by_polygon(baseline, polygon)
    # Expect three segments: outside, inside, outside
    assert len(segments) == 3
    assert segments[1][1] is True
    # intersection segment length approx 222 km
    from hsa_analysis.analysis import geod
    inside_length = geod.geometry_length(LineString(segments[1][0]))
    assert math.isclose(inside_length / 1000, 222, rel_tol=0.05)
