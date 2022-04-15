#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 15:21:10 2021

@author: dabrahamsson
"""

import pandas as pd
import numpy as np

def read_data():
    df = pd.read_csv('expfr_ts_sample3.csv')
    inchi = df.loc[:, 'INCHIKEY']
    dfy = df.loc[:, 'Al_COO':'phenol_noOrthoHbond']
    dfy = np.round(dfy, 0)
    dfy = np.abs(dfy)
    dfy = dfy.astype(int)
    dfy.columns = dfy.columns + '_y'
    df = pd.concat([inchi, dfy], axis=1)
    return(df)

df = read_data()


#uncertainties
def uncertainty(df):
    df = df.set_index('INCHIKEY')
    dft = df.T
    dft = dft.reset_index()
    dft.columns = dft.columns.str.replace('index', 'Fragment')
    dft['Fragment'] = dft['Fragment'].str.replace('_y', '')
    dft = dft.sort_values(by='Fragment')
    dft = dft.set_index('Fragment')
    dft = dft.astype(float)
    
    dfu = pd.read_csv('R2valuesFragmentsTesting.csv')
    dfu.columns = dfu.columns.str.replace('Unnamed: 0', 'Fragment')
    dfu = dfu.sort_values(by='Fragment')
    dfu = dfu.set_index('Fragment')
    dfu = dfu.astype(float)
    
    dfm = pd.read_csv('MAEvaluesFragmentsTesting.csv')
    dfm.columns = dfm.columns.str.replace('Unnamed: 0', 'Fragment')
    dfm = dfm.sort_values(by='Fragment')
    dfm = dfm.set_index('Fragment')
    dfm = dfm.astype(float)
    
    dfut = pd.merge(dft, dfu, on='Fragment')
    dfut = dfut.multiply(dfut["r2"], axis="index")
    dfut = dfut.drop('r2', axis=1)
    
    dfmt = pd.merge(dft, dfm, on='Fragment')
    dfmt = dfmt.multiply(dfmt["mae"]/10, axis="index")
    dfmt = dfmt.drop('mae', axis=1)
    
    #dfg = dfut-dfmt
    dfg = dfmt
    
    dfg = dfg.T
    dfg = dfg.astype(float)
    dfg['error'] = dfg.sum(axis=1)
    dfg = dfg.replace(0, np.NaN)
    dfg['error'] = dfg['error'].fillna(0)
    #dfg['error2'] = dfg['error'] #devide by a factor or 3 to minimize the impact on r2
    
    dfg = dfg.fillna(0)
    dfg = dfg.reset_index()
    dfg.columns = dfg.columns.str.replace('index', 'inchikey_sample')
    dfg = dfg[['inchikey_sample', 'error']]
    return(dfg)

dfg = uncertainty(df)

#read the database 
def read_db():
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

db1, db2, db3 = read_db()

#get the masses for the chemicals in the sample from the database
def get_masses():
    dfs = pd.merge(df, db1, on='INCHIKEY', how='inner')
    
    dfs_inchi = dfs.loc[:, 'INCHIKEY']
    dfs_mass = dfs.loc[:, 'AVERAGE_MASS']
    dfs_bits = dfs.loc[:, 'Al_COO_y':'phenol_noOrthoHbond_y']
    
    dfs1 = pd.concat([dfs_mass, dfs_bits], axis=1)
    dfs2 = pd.concat([dfs_inchi, dfs_mass, dfs_bits], axis=1)
    dfs3 = pd.concat([dfs_inchi, dfs_bits], axis=1)
    return(dfs1, dfs2, dfs3)

dfs1, dfs2, dfs3 = get_masses()

#shorten the database search by searching only for the masses that are present in the sample
def db_search_step1():
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
    
    ds = ds.fillna(0)
    return(dbr1, dbr2, ds)

dbr1, dbr2, ds = db_search_step1()

def db_search_step2(ds):
    colnames = ds.columns.values
    ds = ds.reset_index()
    dff = pd.melt(ds, id_vars=['INCHIKEY'], value_vars=colnames)

    #dff = dff.reset_index()
    dff.columns = dff.columns.str.replace('variable', 'inchikey_db')
    dff.columns = dff.columns.str.replace('INCHIKEY', 'inchikey_sample')
    dff.columns = dff.columns.str.replace('value', 'r')
    return(dff)

dff = db_search_step2(ds)

#filter out mass mismatches
def filter_mismatch():
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

dff1 = filter_mismatch()

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


def show_top_matches(i):
    dft = dg.groupby('inchikey_sample').head(i)
    dft = dft.reset_index(drop=True)
    return(dft)

for i in range(1, 6):
    dft = show_top_matches(i)
    
    # Check for correct matches
    def check_matches(dft):
        dft['match'] = np.where(dft['inchikey_sample'] == dft['inchikey_db'], 1, 0)
        dft = pd.merge(dft, dfg, on='inchikey_sample')
        #Score calculations
        dft['r2'] = dft['r']**2
        dft['Score'] = dft['r2']-dft['error']
        return(dft)
    
    dft = check_matches(dft)

    dft.to_csv('expsimul_ts_sample_3.{}.csv'.format(i))
    print(dft['match'].sum())
    
    
    




























