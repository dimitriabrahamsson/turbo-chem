#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 15:21:10 2021

@author: dabrahamsson
"""

# This file simulates the database search for chemicals with in silico generated 
# physicochemical fingerprints downloaded from the UFZ-LSER database

import pandas as pd
import numpy as np

# Create an in silico sample with 100 chemicals 
def create_sample():
    df1 = pd.read_csv('training_results_exposome_i4.0.csv')
    df2 = pd.read_csv('training_results_exposome_i4.1.csv')
    df3 = pd.read_csv('training_results_exposome_i4.2.csv')
    df4 = pd.read_csv('training_results_exposome_i4.3.csv')
    df5 = pd.read_csv('training_results_exposome_i4.4.csv')
    df = pd.concat([df1, df2, df3, df4, df5], axis=0)
    inchi = df.loc[:, 'INCHIKEY']
    dfy = df.loc[:, 'Al_COO_y':'phenol_noOrthoHbond_y']
    dfy = np.round(dfy, 0)
    dfy = np.abs(dfy)
    dfy = dfy.astype(int)
    df = pd.concat([inchi, dfy], axis=1)
    df = df.sample(n=100)     #Sample 100 random chemicals (rows)
    return(df)

df = create_sample()

# Read the database 
def read_database():
    db = pd.read_csv('BloodExposomeLSERfragment_atoms_R.csv')
    db_inchi = db.loc[:, "INCHIKEY"]
    db_bits = db.loc[:, 'Al_COO':'phenol_noOrthoHbond']
    db_bits.columns = db_bits.columns + '_x'
    db_mass = db.loc[:, 'AVERAGE_MASS']
    df_form = db.loc[:, 'MOLECULAR_FORMULA']
    db1 = pd.concat([db_inchi, db_mass], axis=1)
    db2 = pd.concat([db_inchi, db_bits, db_mass], axis=1)
    db3 = pd.concat([db_inchi, df_form], axis=1)
    return(db1, db2, db3)

db1, db2, db3 = read_database()

# Get the masses for the chemicals in the sample from the database
def get_masses():
    dfs = pd.merge(df, db1, on='INCHIKEY', how='inner')
    dfs_inchi = dfs.loc[:, 'INCHIKEY']
    dfs_mass = dfs.loc[:, 'AVERAGE_MASS']
    dfs_bits = dfs.loc[:, 'Al_COO_y':'phenol_noOrthoHbond_y']
    dfs1 = pd.concat([dfs_mass, dfs_bits], axis=1)
    dfs2 = pd.concat([dfs_inchi, dfs_mass, dfs_bits], axis=1)
    dfs3 = pd.concat([dfs_inchi, dfs_bits], axis=1)
    return (dfs1, dfs2, dfs3)

dfs1, dfs2, dfs3 = get_masses()

def create_matrix():
    #shorten the database search by searching only for the masses that are present in the sample
    dbr1 = db2[db2['AVERAGE_MASS'].isin(dfs1['AVERAGE_MASS'])]
    dbr2 = dbr1.drop('AVERAGE_MASS', axis=1)
    #prepare datasets for correlation matrix
    dfs3.columns = dfs3.columns.str.replace('_y', '')
    dbr2.columns = dbr2.columns.str.replace('_x', '')
    dfs3['INCHIKEY'] = dfs3['INCHIKEY'].astype(str) + '_y'
    dbr2['INCHIKEY'] = dbr2['INCHIKEY'].astype(str) + '_x'
    ds = pd.concat([dfs3, dbr2], axis=0)
    ds = ds.set_index('INCHIKEY')
    ds = ds.T
    ds = ds.corr()
    ds = ds.reset_index()
    ds.columns.name = None
    ds = ds[ds['INCHIKEY'].str.contains('_y')]
    ds = ds.set_index('INCHIKEY')
    ds = ds.loc[:, ds.columns.str.contains('_x')]
    return (ds, dbr1)

ds, dbr1 = create_matrix()

print(ds)

def filterout_mismatches():
    colnames = ds.columns.values
    dsr = ds.reset_index()
    dff = pd.melt(dsr, id_vars=['INCHIKEY'], value_vars=colnames)
    dff.columns = dff.columns.str.replace('variable', 'inchikey_db')
    dff.columns = dff.columns.str.replace('INCHIKEY', 'inchikey_sample')
    dff.columns = dff.columns.str.replace('value', 'r')
    #filter out mass mismatches
    dbrm = dbr1.loc[:, ['INCHIKEY', 'AVERAGE_MASS']]
    dbrm.columns = dbrm.columns.str.replace('INCHIKEY', 'inchikey_db')
    dff['inchikey_db'] = dff['inchikey_db'].str.replace('_x', '')
    dff1 = pd.merge(dff, dbrm, on='inchikey_db')
    dff1.columns = dff1.columns.str.replace('AVERAGE_MASS', 'mass_db')
    dfsm = dfs2.loc[:, ['INCHIKEY', 'AVERAGE_MASS']]
    dfsm.columns = dfsm.columns.str.replace('INCHIKEY', 'inchikey_sample')
    dff1['inchikey_sample'] = dff1['inchikey_sample'].str.replace('_y', '')
    dff1 = pd.merge(dff1, dfsm, on='inchikey_sample')
    dff1.columns = dff1.columns.str.replace('AVERAGE_MASS', 'mass_sample')
    dff1 = dff1[dff1['mass_db'] == dff1['mass_sample']]
    dff1 = dff1.sort_values(by=['inchikey_sample', 'r'], ascending=False)
    return(dff1)

dff1 = filterout_mismatches()

dff1.columns

# Filter our formula mismatches
## This step assumes that the formulas have been assigned correctly
## and that you were able to observe the distinct isotopic patterns of chlorine and bromine
 
def get_formulas():
    dg = pd.merge(dff1, db3, left_on='inchikey_sample', right_on='INCHIKEY', how='left')
    dg.columns = dg.columns.str.replace('MOLECULAR_FORMULA', 'formula_sample')
    dg = dg.drop('INCHIKEY', axis=1)
    dg = pd.merge(dg, db3, left_on='inchikey_db', right_on='INCHIKEY', how='left')
    dg.columns = dg.columns.str.replace('MOLECULAR_FORMULA', 'formula_db')
    dg = dg.drop('INCHIKEY', axis=1)
    dg['formula_agr'] = np.where(dg['formula_sample'] == dg['formula_db'], 1, 0)
    dg = dg[dg['formula_agr'] == 1]
    return(dg)

dg = get_formulas()

# Print 5 files 
# Each time matching to a different number of isomers in the database
# Select how many matches per compound you want to see in range (x, y)

def show_top_matches(i):
    dft = dg.groupby('inchikey_sample').head(i)
    dft = dft.reset_index(drop=True)
    return(dft)

for i in range(1, 6):
    dft = show_top_matches(i)
    
    # Check for correct matches
    def check_matches(dft):
        dft['match'] = np.where(dft['inchikey_sample'] == dft['inchikey_db'], 1, 0)
        return(dft)
    
    dft = check_matches(dft)
    
    dft.to_csv('sim_tr_sample_3.{}.csv'.format(i))
    print(dft['match'].sum())

































