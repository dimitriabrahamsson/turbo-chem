#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 20:15:39 2022

@author: dabrahamsson
"""

import os, glob
import pandas as pd
import shutil, glob

source = '/Users/dabrahamsson/Dropbox/UCSF postdoc/TurboChem experiments/Orbitrap'

destination = '/Users/dabrahamsson/Dropbox/UCSF postdoc/TurboChem experiments/Orbitrap/neg-imputations'
os.chdir(destination)


sol = ['1-octanol', 'butyl acetate', 'chloroform', 'cyclohexane',
       'dichloromethane', 'triolein', 'n-hexane', 'n-octane', 'oleylalcohol',
       'toluene', 'n-undecane']

def combine(sol):
    extension = 'csv'
    all_filenames = [i for i in glob.glob('*.{}'.format(extension))]
    combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames ])
    combined_csv = combined_csv.groupby(['chem_id']).mean()
    combined_csv = combined_csv.drop('Unnamed: 0', axis=1) 
    combined_csv.columns = combined_csv.columns.str.replace('predicted_', '')
    combined_csv = combined_csv[sol]
    combined_csv.columns = 'predicted_' + combined_csv.columns
    return(combined_csv)
    
combined_csv = combine(sol)

combined_csv.to_csv( "neg_combined_csv.csv", encoding='utf-8-sig')

