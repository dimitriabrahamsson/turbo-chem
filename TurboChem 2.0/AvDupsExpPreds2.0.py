# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 13:36:36 2023

@author: dimit
"""

import numpy as np
import pandas as pd


df1 = pd.read_csv('PredsExpData1.csv')
df2 = pd.read_csv('PredsExpData2.csv')
df3 = pd.read_csv('PredsExpData3.csv')
df4 = pd.read_csv('PredsExpData4.csv')
df5 = pd.read_csv('PredsExpData5.csv')

'''
def rounds(df):
    num = df._get_numeric_data()
    num[num < 0] = 0
    df = round(df)
    return(df)

df1 = rounds(df1)
df2 = rounds(df2)
df3 = rounds(df3)
df4 = rounds(df4)f
df5 = rounds(df5)
'''

df = (df1 + df2 + df3 + df4 + df5) / 5
df = df.drop('Unnamed: 0', axis=1)

df['NH'] = df['NH0'] + df['NH1'] + df['NH2']
df = df.drop(['NH0', 'NH1', 'NH2'], axis=1)
#df['bicyclic'] = df['bicyclic'] - 5 # empirical correction of bicyclic - model overpredicts low values 
#df['Al_COO'] = df['Al_COO'] - 1 # empirical correction of bicyclic - model overpredicts low values 
#df['COO'] = df['COO'] - 1 # empirical correction of bicyclic - model overpredicts low values 
#df['COO2'] = df['COO2'] - 1 # empirical correction of bicyclic - model overpredicts low values 

df = df.reindex(sorted(df.columns), axis=1) # Sort df columns alphabetically
df.columns

#df = rounds(df)

print(df['benzene'])


exp = pd.read_excel('ExpDataFragment_atoms2.0.xlsx')
exp.columns.values

e1 = exp.loc[:, 'DTXSID':'Correlation_R2_std2']
e2 = exp.loc[:, 'Al_COO':'urea']
e3 = exp.loc[:, 'H':'Br']

e2['NH'] = e2['NH0'] + e2['NH1'] + e2['NH2']
e2 = e2.reindex(sorted(e2.columns), axis=1) # Sort df columns alphabetically
exp = pd.concat([e1, e2, e3], axis=1)

loc = exp.loc[:, 'DTXSID':'MS_Ready_Monoisotopic_Mass']
dfm = pd.concat([loc, df], axis=1)
dfm.to_csv('AvDupsExpPredsN.csv')

dfmRT = dfm.loc[:, 'Al_COO':'phenol_noOrthoHbond']
dfmRT = dfmRT.set_index(dfm['Preferred_Name'])
dfmRT = dfmRT.T
print(dfmRT)
dfmRT.to_csv('AvDupsT.csv')

expR = exp.loc[:, exp.columns.isin(dfm.columns)]
expR.to_csv('ExpDataFragAtRforCompN.csv')

expRT = expR.loc[:, 'Al_COO':'phenol_noOrthoHbond']
expRT = expRT.set_index(dfm['Preferred_Name'])
expRT = expRT.T
print(expRT)
expRT.to_csv('ExpDataFragAtCompT.csv')
