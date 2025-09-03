# HSAF

Utilities for analyzing maritime baselines and navigation data using the WGS‑84
reference ellipsoid. The code demonstrates how to compute the portion of a
baseline that lies within a specified polygonal region and can serve as a
starting point for more advanced geospatial analysis tools.

## Features
- Calculate length of baseline segments inside a polygon (e.g., territorial sea).
- Split baseline into segments classified as inside or outside the region.
- Simple command line interface for ad‑hoc analysis.

## Installation
```bash
pip install .
```

## Usage
```bash
python -m hsa_analysis.cli --baseline 110 10 112 10 112 12 --polygon 111 9 113 9 113 11 111 11
```

## Testing
```bash
pytest
```

The repository retains a historical MATLAB script (`Hankel_Filter.m`) that is
unrelated to the new baseline analysis utilities but preserved for reference.
