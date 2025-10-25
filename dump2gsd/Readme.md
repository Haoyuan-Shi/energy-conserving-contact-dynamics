# dump2gsd: Convert LAMMPS Dump to GSD for OVITO Visualization

This script converts LAMMPS custom dump files (with quaternion orientation data) into `.gsd` trajectory files for visualization in [OVITO](https://www.ovito.org/). It is tailored for simulations of rigid cube particles represented using body-based pair styles.


## Features

- Parses LAMMPS dump files with quaternion format
- Wraps particle positions with periodic boundary conditions (PBC)
- Constructs cube-shaped particles with proper vertex/face geometry
- Creates a HOOMD `.gsd` trajectory compatible with OVITO


## Requirements

- Python 3
- `hoomd` (tested with HOOMD-blue v4+)
- `gsd`
- `numpy`
- `pandas`

Install via pip if needed:
```python
pip install hoomd gsd numpy pandas
```

## Input Dump Format

The script expects a custom LAMMPS dump file generated using the following commands:
```lammps
compute 10 all property/atom x y z quatw quati quatj quatk
dump mydump all custom 100 dump.quat id c_10[1] c_10[2] c_10[3] c_10[4] c_10[5] c_10[6] c_10[7]
```

## Usage

Run the script with the dump file as input:
```bash
python dump2gsd.py dump.quat
```
This will generate:

- `trajectory.gsd` — a HOOMD-compatible trajectory file viewable in OVITO

Intermediate files like `cubes.gsd` will be automatically removed after use.


## Cube Geometry

The cube particles are defined using vertex and face lists. Each cube has edge length 1.0 and a roundness radius of 0.15.


## Output

`trajectory.gsd`: Contains per-frame cube positions, orientations, and box dimensions.

Can be directly opened with OVITO or used in further HOOMD analysis.


## Notes

- The script assumes particles are labeled with the type `cube`.
- Default velocity is set to zero for all particles.
- Only positions and quaternions are read from the dump file; no forces or other properties are preserved.


## License

MIT License


## Author

Haoyuan Shi
