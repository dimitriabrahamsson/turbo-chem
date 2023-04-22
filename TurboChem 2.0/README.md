## Alright, welcome to TurboChem 2.0!

### This tutorial will guide you through the steps of using TurboChem and reproducing the data in the study "Extracting structural information from physicochemical properties - A new approach in structure elucidation for non-targeted analysis. 

#### To keep things simple, let's start with the negative ionization data (ESI-) from the QTOF method.
The script you want to start with is called DataProcNeg_step1.py. This script will read two files: 1) NegAreasExp1.xlsx and 2) ENTACT_504.xlsx. The first file contains the peak areas for all the chemical features from experiment 1 and the second file contains the chemical names and monoisotopic masses for the compounds contained in the ENTACT mixture 504. Run the script to clean the data and prepare them for the next step.

#### The next script is DataProcNeg_step2.py. 
This script will read the file you generated from the previous script and the database TurboChemDB5.2.csv to impute the missing values based on the process described in the paper.

