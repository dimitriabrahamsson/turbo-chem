## Alright, welcome to TurboChem 2.0!

### This tutorial will guide you through the steps of using TurboChem and reproducing the data in the study "Extracting structural information from physicochemical properties - A new approach in structure elucidation for non-targeted analysis. 

#### To keep things simple, let's start with the negative ionization data (ESI-) from the QTOF method.
The script you want to start with is called DataProcNeg_step1.py. This script will read two files: 1) NegAreasExp1.xlsx and 2) ENTACT_504.xlsx. The first file contains the peak areas for all the chemical features from experiment 1 and the second file contains the chemical names and monoisotopic masses for the compounds contained in the ENTACT mixture 504. Run the script to clean the data and prepare them for the next step.

#### The next script is DataProcNeg_step2.py. 
This script will read the file you generated from the previous script and the database TurboChemDB5.2.csv to impute the missing values based on the process described in the paper.

#### Now let's combine all the imputation files into one
Move all generated files from the previous step into one folder. Run script DataProcNeg_step3.py inside that folder to combine all files into one by averaging the imputations.

#### The last step in data processing
In this final step run script DataProcNeg_step4.py to read your original peak areas and fill in the gaps using the imputations.

#### Now let's make some structure predictions
First we need to run the script called ANN_structure_pred_expData3.1.py. This script contains the artificial neural net that we will train to make predictions. The first part of the script reads the database file TurboChemDB5.2.csv and prepares the training and testing sets. Then it reads the experimental data from a file called ExpDataFragment_atoms2.0.xlsx, which contains the peak areas and the Ksw from the experiments. 

Once the model is trained and it will make predictions for the training, testing and experimental sets. The output is shown in the files PredsExpData5.csv, PredsTrain5.csv, PredsTest5.csv. 

You might want to repeat the process 5 times and take the average of the 5 predictions. In this case run script AvDupsExpPreds2.0.py to average the predictions. 