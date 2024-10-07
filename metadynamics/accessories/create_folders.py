import os
import shutil

# Define sigma values
temp_values = [50, 100, 150, 200, 250, 300, 350, 400, 450]

# Load the content of the files

with open('input.xml', 'r') as file:
    input_xml_content = file.read()

with open('run.sh', 'r') as file:
    run_sh_content = file.read()

with open('run-mace.py', 'r') as file:
    mace_content = file.read()
with open('plumed-md.dat', 'r') as file:
    plumed_content = file.read()

# Folder to save individual runs
base_folder = 'individual_dist'
if not os.path.exists(base_folder):
    os.mkdir(base_folder)

# Initialize the content for run_all.sh
run_all_content = "#!/bin/bash\n\n"

for sigma in temp_values:
    # Create a unique folder for each sigma value
    folder_name = f'{base_folder}/{sigma}K'
    os.makedirs(folder_name, exist_ok=True)

    # Modify the plumed-md.dat with the current sigma value

    
    # Modify the input.xml with the updated driver name (based on sigma)
    driver_name = f'driver-works-{sigma}'
    modified_input_xml = input_xml_content.replace('driver-combo', driver_name)
    modified_input_xml = modified_input_xml.replace('400', str(sigma))
    modified_plumed = plumed_content.replace('400', str(sigma))
    # Modify run.sh with the updated driver name (based on sigma)
    modified_run_sh = run_sh_content.replace('driver-combo', driver_name)

    # Modify run-mace.py with the updated driver name (based on sigma)
    modified_mace = mace_content.replace('driver-combo', driver_name)

    # Save the modified files in the new folder

    with open(f'{folder_name}/input.xml', 'w') as file:
        file.write(modified_input_xml)

    with open(f'{folder_name}/run.sh', 'w') as file:
        file.write(modified_run_sh)

    with open(f'{folder_name}/run-mace.py', 'w') as file:
        file.write(modified_mace)
    with open(f'{folder_name}/plumed-md.dat', 'w') as file:
        file.write(modified_plumed)


    # Copy other necessary files like init.xyz to the folder
    shutil.copy('init.xyz', folder_name)
    # Add the command to run the simulation in parallel in the new folder
    run_all_content += f"(cd {folder_name} && bash run.sh) &\n"
    
    print(f'Prepared folder for sigma {sigma}')

# Write the run_all.sh script
run_all_content += "wait\n"
with open(f'{base_folder}/run_all.sh', 'w') as file:
	file.write(run_all_content)

# Make run_all.sh executable
os.chmod(f'{base_folder}/run_all.sh', 0o755)

print("run_all.sh script generated successfully!")
