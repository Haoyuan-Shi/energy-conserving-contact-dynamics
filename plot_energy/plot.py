import numpy as np
import matplotlib.pyplot as plt
import sys

plt.rc('font', size=8)
plt.rc('lines', linewidth=1)
plt.rc('axes', linewidth=1)
plt.rc('font', family='serif')
plt.rcParams.update({'figure.autolayout': True})

# Check for correct usage
if len(sys.argv) != 2:
    print("Usage: python plot.py <number_of_atoms>")
    sys.exit(1)

# Read number of atoms from command line
try:
    n = int(sys.argv[1])
except ValueError:
    print("Error: number_of_atoms must be an integer.")
    sys.exit(1)

print(f"Number of atoms: {n}")

eng = np.loadtxt("./energy.dat")
a1 = (3*n-3)/2
a2 = (6*n-3)/2
fig, ax = plt.subplots(figsize=(3.2, 3.45))  # Use plt.subplots() instead of plt.Figure()

steps = np.arange(0, len(eng[:,0])) * 100

total_energy = eng[:, 3]+eng[:, 5]*a2
ax.plot(steps, eng[:, 5]*a2 - eng[:, 1]*a1, color='green', linestyle=':', label='Rot. kinetic energy')
ax.plot(steps, eng[:, 1]*a1, color='purple', linestyle=':', label='Trans. kinetic energy')
ax.plot(steps, eng[:, 3], color='firebrick', linestyle='--', label='Potential energy')
ax.plot(steps, eng[:, 5]*a2, color='midnightblue', linestyle='--', label='Kinetic energy')
ax.plot(steps, total_energy, color='black', linestyle='-', label='Total energy')

ax.set_xlabel('Step')
ax.set_ylabel('Energy')
# ax.set_ylim([-70,700])
ax.legend(bbox_to_anchor=(-0.2, 1.0), loc='lower left', ncols=2, frameon=False)

fig.savefig('energy-conservation.png', bbox_inches='tight', dpi=1000)
print(total_energy[-1] - total_energy[0])
