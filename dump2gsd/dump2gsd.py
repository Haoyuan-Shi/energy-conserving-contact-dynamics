import os
import sys
import numpy as np
import hoomd
import gsd.hoomd
import pandas as pd

vertices = [
    (0.5, 0.5, 0.5),
    (0.5, -0.5, 0.5),
    (-0.5, -0.5, 0.5),
    (-0.5, 0.5, 0.5),
    (0.5, 0.5, -0.5),
    (0.5, -0.5, -0.5),
    (-0.5, -0.5, -0.5),
    (-0.5, 0.5, -0.5),
]
faces = [
    (2, 1, 0),
    (0, 3, 2),
    (5, 6, 4),
    (4, 6, 7),
    (1, 2, 5),
    (5, 2, 6),
    (4, 7, 3),
    (3, 0, 4),
    (0, 1, 5),
    (5, 4, 0),
    (3, 7, 6),
    (6, 2, 3)
]

radius = 0.15

def read_lammps_custom_quat(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    frames = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            timestep = int(lines[i+1].strip())
            num_atoms = int(lines[i+3].strip())
            box_bounds = [list(map(float, lines[i+5].strip().split())),
                          list(map(float, lines[i+6].strip().split())),
                          list(map(float, lines[i+7].strip().split()))]
            atom_header = lines[i+8].strip().split()[2:]
            atom_lines = lines[i+9:i+9+num_atoms]
            
            # Parse atom data
            data = [list(map(float, line.strip().split())) for line in atom_lines]
            df = pd.DataFrame(data, columns=atom_header)
            
            # Store all data for the frame
            frames.append({
                "timestep": timestep,
                "box_bounds": box_bounds,
                "atoms": df
            })
            
            i = i + 9 + num_atoms
        else:
            i += 1

    return frames

# Load the dump file
if len(sys.argv) != 2:
    print("Usage: python dump2gsd.py <LAMMPS CUSTOM DUMP FILE id x y z quatw quati quatj quatk>")
    sys.exit(1)

filename = sys.argv[1]
frames = read_lammps_custom_quat(filename)

def wrap_position(coord, lo, hi):
    length = hi - lo
    return lo + ((coord - lo) % length)

timestep = []
box_bounds = []
atom = []
for frame in frames:
    timestep.append(frame["timestep"])
    box_bounds.append(frame["box_bounds"])
    atoms = frame["atoms"]
    atoms = atoms.rename(columns={
        'c_10[1]': 'x',
        'c_10[2]': 'y',
        'c_10[3]': 'z',
        'c_10[4]': 'quatw',
        'c_10[5]': 'quati',
        'c_10[6]': 'quatj',
        'c_10[7]': 'quatk',
    })
    atoms["id"] = atoms["id"].astype(int)

    atoms["x"] = wrap_position(atoms["x"], frame["box_bounds"][0][0], frame["box_bounds"][0][1])
    atoms["y"] = wrap_position(atoms["y"], frame["box_bounds"][1][0], frame["box_bounds"][1][1])
    atoms["z"] = wrap_position(atoms["z"], frame["box_bounds"][2][0], frame["box_bounds"][2][1])

    atom.append(atoms)

n_particles = len(atom[0])
particle_ids = atom[0]["id"]

if os.path.exists('trajectory.gsd'):
    os.remove('trajectory.gsd')

for ii in range(len(frames)):
    cell = hoomd.md.nlist.Cell(buffer=0.1)
    alj = hoomd.md.pair.aniso.ALJ(nlist=cell, default_r_cut=0)
    alj.params[('cube', 'cube')] = dict(
        epsilon=0, 
        sigma_i=2*radius, 
        sigma_j=2*radius, 
        alpha=0,
        contact_ratio_i=1,
        contact_ratio_j=1
    )
    
    alj.shape['cube'] = dict(
        vertices=vertices,
        faces=faces
    )
    integrator = hoomd.md.Integrator(dt=1.0E-8, integrate_rotational_dof=True)
    integrator.forces.append(alj)    
        
    cpu = hoomd.device.CPU()
    simulation = hoomd.Simulation(device=cpu, seed=0)
    simulation.operations.integrator = integrator
    
    position = []
    orientation = []
    for pid in particle_ids:
        row = atom[ii][atom[ii]["id"] == pid]
        position.append((row["x"], row["y"], row["z"]))
        orientation.append((row["quatw"], row["quati"], row["quatj"], row["quatk"]))
    
    frame = gsd.hoomd.Frame()
    frame.particles.N = n_particles
    frame.particles.position = position
    frame.particles.orientation = orientation
    frame.particles.typeid = [0] * n_particles
    frame.particles.types = ['cube']
    frame.configuration.box = [box_bounds[ii][0][1]-box_bounds[ii][0][0], box_bounds[ii][1][1]-box_bounds[ii][1][0], box_bounds[ii][2][1]-box_bounds[ii][2][0], 0, 0, 0]

    velocity = []
    for pid in particle_ids:
        velocity.append((0.0, 0.0, 0.0))
    frame.particles.velocity = velocity
    
    with gsd.hoomd.open(name='cubes.gsd', mode='w') as f:
        f.append(frame)
    
    simulation.create_state_from_gsd(filename='cubes.gsd')
    logger = hoomd.logging.Logger()
    logger.add(
        alj, 
        quantities=[
            'type_shapes'
        ]
    )
    gsd_writer = hoomd.write.GSD(
        filename='trajectory.gsd', 
        logger=logger,
        trigger=hoomd.trigger.Periodic(1), 
        mode='ab', 
        filter=hoomd.filter.All(), 
    )
    simulation.operations.writers.append(gsd_writer)
    simulation.run(1)
    gsd_writer.flush()

if os.path.exists('cubes.gsd'):
    os.remove('cubes.gsd')