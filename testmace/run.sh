IPI=i-pi
PYTHON=python

${IPI} input.xml &> log.i-pi & 

sleep 10

${PYTHON} run-mace.py & 


wait
