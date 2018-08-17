pip install biopython
pip install pandas
pip install synbiochem-py

export PYTHONPATH=$PYTHONPATH:.

python \
	pcr/pcr_genie.py \
	https://ice.synbiochem.co.uk \
	ICE_USERNAME \
	ICE_PASSWORD \
	ice_ids.txt \
	gaattcaaaagatcttttaagaag \
	ttactcgagtttggatcc \
	out \
	lengths.csv