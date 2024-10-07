import numpy as np

# Load the data from the .dat file, skipping the first line (header)
data = np.loadtxt('fullerenes/new_1/fullerene_3/surface/rotate/energy_landscape.dat', skiprows=1)

# Split data into coordinates and energy values
coords = data[:, :3]
energies = data[:, 3]

# Write to VTK format
with open('fullerenes/new_1/fullerene_3/surface/rotate/energy_grid.vtk', 'w') as f:
    f.write("# vtk DataFile Version 3.0\n")
    f.write("VTK file with energy grid data\n")
    f.write("ASCII\n")
    f.write("DATASET UNSTRUCTURED_GRID\n")
    
    # Write points (x, y, z coordinates)
    num_points = len(coords)
    f.write(f"POINTS {num_points} float\n")
    for coord in coords:
        f.write(f"{coord[0]:.10f} {coord[1]:.10f} {coord[2]:.10f}\n")
    
    # Write cells (each point is a vertex)
    f.write(f"CELLS {num_points} {2 * num_points}\n")
    for i in range(num_points):
        f.write(f"1 {i}\n")
    
    f.write(f"CELL_TYPES {num_points}\n")
    for _ in range(num_points):
        f.write("1\n")  # VTK_VERTEX for all points

    # Write scalar data (energy)
    f.write(f"POINT_DATA {num_points}\n")
    f.write("SCALARS energy float 1\n")
    f.write("LOOKUP_TABLE default\n")
    for energy in energies:
        f.write(f"{energy:.10f}\n")
