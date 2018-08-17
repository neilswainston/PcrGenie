c:\Python27\Scripts\pip install biopython
c:\Python27\Scripts\pip install pandas
c:\Python27\Scripts\pip install synbiochem-py

set PYTHONPATH=%PYTHONPATH%;.

python pcr/pcr_genie.py https://ice.synbiochem.co.uk ICE_USERNAME ICE_PASSWORD ice_ids.txt gaattcaaaagatcttttaagaag ttactcgagtttggatcc out lengths.csv

pause