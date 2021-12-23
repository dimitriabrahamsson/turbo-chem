## turbo-chem
### In silico structure prediction for non-targeted analysis with high-resolution mass spectrometry - From physicochemical properties to molecular structures

#### Preparing the data  
1. Start with the .cvs file called BloodExposomeLSERfinal.csv  
2. Run blood_exposome_frag_atoms.py to generate BloodExposomeLSERfragment_atoms_R.csv  

#### Convert the physicochemical fingerprints to RDKit fragments
Run ANN_structure_pred.py to generate the result files for the training and testing sets training_results_exposome_i1.0.csv and test_results_exposome_i1.0.csv

#### Simulate the database search (for evaluation only)
Run DB_search_sim_fr.py to simulate the database searching using data from the training and testing sets. This step is meant to be used as a preliminary step for evaluating the algorithm

#### Search the database with fragments generated from experimental physicochemical fingerprints
Run DB_search_exp_sim_fr.py to search the database using RDKit fragments that were generated with experimental physicochemical fingerprints. To generate the fragments use the blood_exposome_frag_atoms.py script.


