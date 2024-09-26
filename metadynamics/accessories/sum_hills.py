import os
import subprocess

# List of folders to run the command in
folders = [200,300,400] #this is for different temperatures

# Define the command to be run
command = [
    "plumed", "sum_hills",
    "--hills", "HILLS-md",
    "--min", "0.0",
    "--max", "2.3",
    "--bin", "100",
    "--outfile", "FES-md2",
  "--stride", "100", #can delete this - this is to ouput multiple times at different timesteps to check convergence
	"--mintozero"
]
# Loop over each folder
for folder in folders:
    folder_path = str(str(folder) + "K")  # Convert folder number to string
    if os.path.exists(folder_path):  # Check if the folder exists
        print(f"Running command in folder: {folder_path}")
        # Change to the folder
        os.chdir(folder_path)

        # Run the plumed sum_hills command
        result = subprocess.run(command, capture_output=True, text=True)

        # Print the output and error if any
        print(f"Output in folder {folder_path}:\n{result.stdout}")
        if result.stderr:
            print(f"Error in folder {folder_path}:\n{result.stderr}")

        # Go back to the original directory
        os.chdir("..")
    else:
        print(f"Folder {folder_path} does not exist.")
