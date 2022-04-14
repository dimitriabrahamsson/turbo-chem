## Welcome to turbo-chem!
### Turbo-chem is an in silico structure prediction algorithm designed to propose molecular structures for chemical features detected  in non-targeted analysis with high-resolution mass spectrometry
#### Turbo-chem is based on the idea that different molecular structures have distinctly different physicochemical properties and when these properties are put together they form a unique fingerpring - a "physicochemical fingerprint".

#### The following steps will help you run turbo-chem on your computer and make structure predictions for a set of chemicals.

#### Preparing the data  
1. Start with the .cvs file called BloodExposomeLSERfinal.csv  
2. Run blood_exposome_frag_atoms.py to generate BloodExposomeLSERfragment_atoms_R.csv  

#### Convert the physicochemical fingerprints to RDKit fragments
Run ANN_structure_pred.py to generate the result files for the training and testing sets training_results_exposome_i1.0.csv and test_results_exposome_i1.0.csv

#### Simulate the database search (for evaluation only)
Run DB_search_sim_fr.py to simulate the database searching using data from the training and testing sets. This step is meant to be used as a preliminary step for evaluating the algorithm

#### Search the database with fragments generated from experimental physicochemical fingerprints
Run DB_search_exp_sim_fr.py to search the database using RDKit fragments that were generated with experimental physicochemical fingerprints. To generate the fragments use the blood_exposome_frag_atoms.py script.


