#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 12:26:02 2020

@author: dabrahamsson
"""

import rdkit
from rdkit import Chem
import rdkit.Chem.Fragments
import rdkit.Chem.FunctionalGroups
import rdkit.Chem.Draw
import pandas as pd
import numpy as np

df = pd.read_csv('BloodExposomeLSERfinal.csv')

#df = df.iloc[:10] 

mols = []
for i in df['QSAR_READY_SMILES']:
    mols.append(Chem.AddHs(Chem.MolFromSmiles(i))) 
    
from fragmentation import fragments
from fragmentation import atom_numbers

count = 0
df_fragments = pd.DataFrame()
for mol in mols:
    df_frag = fragments(mol)
    if count % 100 == 0:
        print(count)
    df_fragments = df_fragments.append(df_frag)
    count +=1 

df_fragments = df_fragments.reset_index(drop=True)
df_fragments

df_fragments.to_csv('BloodExposomeLSERfragments.csv')
df_fr = pd.concat([df, df_fragments], axis=1)
df_fr.to_csv('BloodExposomeLSERfr.csv')


count = 0
df_atoms = pd.DataFrame()
for mol in mols:
    df_atom = atom_numbers(mol)
    if count % 100 == 0:
        print(count)
    df_atoms = df_atoms.append(df_atom)
    count +=1 

df_atoms = df_atoms.reset_index(drop=True)
df_atoms

df_atoms.to_csv('BloodExposomeLSERatoms.csv')
df_at = pd.concat([df, df_atoms], axis=1)
df_at.to_csv('BloodExposomeLSERat.csv')

def feature_extration(df):
    df = df.replace(0, np.NaN)
    column_counts = df.count()
    column_counts = column_counts[column_counts>(len(df)*0.10)]
    column_counts.index
    df = df.loc[:, df.columns.isin(column_counts.index)]
    dfR = df.fillna(0)
    return (dfR)

df_fragR = feature_extration(df_fragments)
df_fragR

df_frag_at = pd.concat([df_fragR, df_atoms], axis=1)

df_frag_at.columns.values

df1 = pd.concat([df_frag_at, df], axis=1)

df1.to_csv('BloodExposomeLSERfragment_atoms_R.csv')
