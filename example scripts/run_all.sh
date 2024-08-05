#!/bin/bash

# Declare directories of individual simulations, input driver names, steps, temperature, and beads
declare -A directories=(
    ["C60"]="driver_4 200000 500 1"
    ["C180"]="driver_5 300000 600 2"
)

# Function to run job in a directory
run_job() {
    dir=$1
    driver=$2
    steps=$3
    temp=$4
    beads=$5
    echo "Running job in $dir with driver $driver, steps=$steps, temperature=$temp, nbeads=$beads"
    cd $dir

    # Ensure run.sh is executable
    chmod +x run.sh

    # Create input file with parameters
    cat > input.xml <<EOF
<simulation verbosity='medium'>
  <output prefix='simulation'>
    <properties filename='out' stride='50' flush='10'>  [ step, time{picosecond}, conserved, temperature{kelvin}, potential ] </properties>
    <trajectory filename='pos' stride='50' flush='100' format='ase'> positions </trajectory>
    <checkpoint stride='200'/>
  </output>
  <total_steps> ${steps} </total_steps>
  <prng>
    <seed>32415</seed>
  </prng>
  <ffsocket name='maceoff23' mode='unix' pbc='false'>
    <address> ${driver} </address>
  </ffsocket>
  <system>
    <initialize nbeads='${beads}'>
      <file mode='xyz'> init.xyz </file>
      <velocities mode='thermal' units='kelvin'> 400 </velocities>
    </initialize>
    <forces>
      <force forcefield='maceoff23' weight='1'> </force>
    </forces>
    <motion mode='dynamics'>
      <fixcom> True </fixcom>
      <dynamics mode='nvt'>
        <timestep units='femtosecond'> 0.5</timestep>
        <thermostat mode='svr'>
          <tau units='femtosecond'> 100 </tau>
        </thermostat>
      </dynamics>
    </motion>
    <ensemble>
      <temperature units='kelvin'> ${temp} </temperature>
    </ensemble>
  </system>
</simulation>
EOF

    # Pass the driver name to run.sh
    ./run.sh $driver &

    cd -  # Go back to the original directory
}

# Run jobs in parallel
for dir in "${!directories[@]}"; do
    params=(${directories[$dir]})
    run_job $dir "${params[0]}" "${params[1]}" "${params[2]}" "${params[3]}" &
done

# Wait for all background jobs to complete
wait

echo "All tasks completed."
