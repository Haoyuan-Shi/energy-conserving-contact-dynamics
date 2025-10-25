# Energy-Conserving Contact Dynamics of Nonspherical Rigid-Body Particles

This repository provides a fully enhanced rigid-body particle simulation for LAMMPS (22 Jul 2025). Based on [Langston et al., 2010](https://link.springer.com/article/10.1007/s10035-010-0217-4), this framework resolve discontinuous forces and expand the model to support complex particle interactions and energy-conserving dynamic behaviors.

## 🚀 Contribution

This implementation introduces key improvements for handling rigid body interactions in simulations, focusing on **energy conservation**.

### Energy Conservation via Conservative Forces

To ensure **energy-conserving dynamics**, interactions have been modified to eliminate **discontinuous forces**. Smooth contact transitions are achieved by refining the contact detection logic based on geometric relevance, resulting in stable and physically consistent simulations.

### Contact Points Detection

**2D Interactions**  
🔹 **Vertex → Boundary**

**3D Interactions**  
🔹 **Vertex → Surface**  
🔹 **Edge → Edge**

These refinements improve the physical fidelity of contact resolution in rigid-body simulations, especially for complex polyhedral interactions.

## Visual Overview

<table style="border-collapse: collapse; text-align: center;">
  <tr>
    <th>Model (125 cubes with density: 0.125)</th>
    <th>Trajectory</th>
    <th>Radial Distribution Function</th>
    <th>Velocity Auto-Correlation</th>
  </tr>
  <tr>
    <td><img src="image/energy.png" alt="Energy" width="300"/></td>
    <td><img src="image/cubes.gif" alt="Trajectory" width="180"/></td>
    <td><img src="image/rdf.png" alt="Radial Distribution Function" width="300"/></td>
    <td><img src="image/vacf.png" alt="Velocity Auto-Correlation Function" width="300"/></td>
  </tr>
</table>

## 🔧 Installation

The codes have been tested with [LAMMPS (22 Jul 2025)](https://download.lammps.org/tars/index.html). Compatibility with newer versions is not guaranteed.

#### Steps

1. **Download and extract LAMMPS**  
   Get the source from the [LAMMPS download page](https://download.lammps.org/tars/index.html).

2. **Copy source files**  
   Copy the following files into the `lammps/src/BODY/` directory:
   - `pair_body_rounded_polyhedron.cpp`
   - `pair_body_rounded_polyhedron.h`

3. **Build LAMMPS with the BODY package**  
   Follow the official instructions to:
   - [Build LAMMPS](https://docs.lammps.org/Build.html)
   - [Include the BODY package](https://docs.lammps.org/Build_package.html)

## 🧭 Visualization

The script `~/dump2gsd/dump2gsd.py` can be used to convert a LAMMPS custom dump file into a `.gsd` file, which can be visualized in [OVITO](https://www.ovito.org/).

> **Note:** The dump command in your LAMMPS input script **must match exactly** the following format:

```lammps
compute 10 all property/atom x y z quatw quati quatj quatk
dump mydump all custom 100 dump.quat id c_10[1] c_10[2] c_10[3] c_10[4] c_10[5] c_10[6] c_10[7]
```

## 📈 Plot Energy

To generate energy plots from your LAMMPS simulation:

1. Copy the contents of the `~/plot_energy` directory into the folder containing your `log.lammps` file.
2. Run the script using:

  ```bash
  bash energy.sh
  ```

## 🧪 Test
A sample system containing 125 cubes is available in the `~/test` directory.
