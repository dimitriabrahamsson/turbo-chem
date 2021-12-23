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

# Creating an in silico sample with 100 chemicals 
def create_sample():
    df = pd.read_csv('test_results_exposome_i1.0.csv')
    inchi = df.loc[:, 'INCHIKEY']
    dfy = df.loc[:, 'Al_COO_y':'phenol_noOrthoHbond_y']
    dfy = np.round(dfy, 0)
    dfy = np.abs(dfy)
    dfy = dfy.astype(int)
    df = pd.concat([inchi, dfy], axis=1)
    df = df.sample(n=100)     #Sample 100 random chemicals (rows)
    return (df)

df = create_sample()

# Read the database 
def read_database():
    db = pd.read_csv('BloodExposomeLSERfragment_atoms_R.csv')
    db_inchi = db.loc[:, "INCHIKEY"]
    db_bits = db.loc[:, 'Al_COO':'phenol_noOrthoHbond']
    db_bits.columns = db_bits.columns + '_x'
    db_mass = db.loc[:, 'AVERAGE_MASS']
    db1 = pd.concat([db_inchi, db_mass], axis=1)
    db2 = pd.concat([db_inchi, db_bits, db_mass], axis=1)
    return(db1, db2)

db1, db2 = read_database()

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

def filterout_missmatches():
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

dff1 = filterout_missmatches()

# Select how many matches per compound you want to see .head(x)
def show_top_matches():
    dft = dff1.groupby('inchikey_sample').head(1)
    dft = dft.reset_index(drop=True)
    return(dft)

dft = show_top_matches()


# Check for correct matches
def check_matches(dft):
    dft['inchikey_sample'] = dft['inchikey_sample'].str.replace('_y', '')
    dft['inchikey_db'] = dft['inchikey_db'].str.replace('_x', '')
    dft['match'] = np.where(dft['inchikey_sample'] == dft['inchikey_db'], 1, 0)
    print(dft)
    print(dft['match'].sum())
    return(dft)

dft = check_matches(dft)

dft.to_csv('ts_sample 1.0.csv')
print(dft['match'].sum())






























