import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

def distance_vs_time(simulation_path, fullerene_path):
    # Load your main trajectory and topology files
    u = mda.Universe(simulation_path, format='XYZ')

    # Load the fullerene XYZ file
    fullerene_u = mda.Universe(fullerene_path)

    # Select all atoms in the fullerene
    fullerene = fullerene_u.select_atoms('name C')

    # Select the oxygen atom of the water molecule
    water_oxygen = u.select_atoms('name O')
    water_oxygen = water_oxygen[[-1]][0]

    def calculate_com(atoms):
        masses = atoms.masses
        positions = atoms.positions
        com = np.average(positions, weights=masses, axis=0)
        return com

    distances = []
    times = []
    num_frames = len(u.trajectory)
    com_fullerene = calculate_com(fullerene)
    
    print(f"Number of frames in trajectory: {num_frames}")
    for ts in u.trajectory:

        pos_water_oxygen = water_oxygen.position

        # Compute distance between COM of fullerene and water oxygen atom
        distance = np.linalg.norm(com_fullerene - pos_water_oxygen)
        distances.append(distance)
        times.append(ts.time)

    # Convert lists to numpy arrays
    distances = np.array(distances)
    times = np.array(times)
    return distances, times
    # Plot distances vs time
def plot_distance_vs_time(simulation_path, fullerene_path, label):
    distances, times = distance_vs_time(simulation_path, fullerene_path)
    plt.plot(times, distances, marker='o', linestyle='-', label=label)
    plt.xlabel('Time (ps)')
    plt.ylabel('Distance (Å)')
    plt.title('Distance of Water Molecule from COM of Fullerene vs Time')
def plot_rdf(simulation_path, fullerene_path, label, bin_width=0.1, r_min = 0.0, r_max=4, sigma=2.0):
    u = mda.Universe(simulation_path, format='XYZ')
    distances, _ = distance_vs_time(simulation_path, fullerene_path)

    n_bins = int(r_max / bin_width)
    bin_edges = np.linspace(r_min, r_max, n_bins + 1)
    rdf, _ = np.histogram(distances, bins=bin_edges, density=True)

    # Normalize RDF
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    volume = (4/3) * np.pi * (bin_edges[1:]**3 - bin_edges[:-1]**3)
    rdf /= volume

    # Apply Gaussian smoothing
    rdf_smooth = gaussian_filter1d(rdf, sigma=sigma)

    # Plot RDF
    plt.plot(bin_centers, rdf_smooth, label=label)
    plt.xlim(r_min, r_max)
    plt.xlabel('Distance (Å)')
    plt.ylabel('g(r)')
    plt.title('RDF of a water molecule from the COM of a fullerene')
    

# Call the function with different simulation paths
for a in range (0, 7):
    plot_rdf(rf"C:\Users\akrya\UROP2024\fullerenes\new_1\mCE_OFF\8bead\sim.pos_{a}.extxyz", r"C:\Users\akrya\UROP2024\fullerenes\new_1\fullerene_only.xyz", f'bead {a}')
for a in ['400K']:
    plot_rdf(rf"C:\Users\akrya\UROP2024\fullerenes\new_1\mCE_OFF\{a}\simulation.pos_0.extxyz", r"C:\Users\akrya\UROP2024\fullerenes\new_1\fullerene_only.xyz", f'{a}')
#plot_rdf(rf"C:\Users\akrya\UROP2024\fullerenes\new_1\mace_mp\simulation.pos_0.extxyz", r"C:\Users\akrya\UROP2024\fullerenes\new_1\fullerene_only.xyz",'MACE-off clasical')
#plot_distance_vs_time('fullerenes/new_1/mce_off/simulation.pos_0.extxyz', 'fullerenes/new_1/fullerene_only.xyz', 'Simulation 1')
plt.legend()
plt.show()
#plot_rdf('fullerenes/new_1/mce_off/simulation.pos_0.extxyz', 'fullerenes/new_1/fullerene_only.xyz', 'Simulation 1:Mace-off')
#plot_rdf('fullerenes/new_1/simulation.pos_0.extxyz', 'fullerenes/new_1/fullerene_only.xyz', 'Simulation 2: Mace-MP')
#plt.legend()
#plt.show()
