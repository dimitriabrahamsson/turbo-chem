#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 09:48:53 2022

@author: dabrahamsson
"""

import numpy as np
import pandas as pd
import seaborn as sns

def readData():
    df = pd.read_excel('Area_0_20227251754_pos.xlsx')
    entact = pd.read_excel('ENTACT_504.xlsx') 
    df1 = df.loc[:, 'Average Rt(min)':'Average Mz']
    df2 = df.loc[:, 'pos_DIWB_1':'pos_no injection_4']
    df3 = df.loc[:, 'pos_1-octl_1':'pos_undecW_2']
    df4 = df.loc[:, 'Average':'Stdev']
    return(entact, df1, df2, df3, df4)

entact, df1, df2, df3, df4 = readData()

def massCorrect(df):
    # Correct MS1 m/z to account for M+H
    df['Average Mz'] = df['Average Mz'] - 1.0078250319
    return(df)
    
df1 = massCorrect(df1)

def averages(df):
    # Take the averages for every two samples
    column_names = df.columns.values
    column_names = pd.DataFrame({'sample_names': column_names})
    column_names.sample_names = column_names.sample_names.str.replace('_1', '')
    column_names.sample_names = column_names.sample_names.str.replace('_2', '')
    column_names = column_names.sample_names.drop_duplicates()
    column_names = column_names.reset_index(drop=True)
    
    df = df.groupby(np.arange(len(df.columns))//2, axis=1).mean()
    df.columns = column_names
    return(df)

df3 = averages(df3)

def missingArea(df):
    loc1 = df[df.columns[::2]]
    loc2 = df.loc[:, df.columns.str.contains('W')]
    loc2.columns = loc2.columns.str.replace('W', '')
    loc3 = loc2-loc1
    loc3 = loc3.add_suffix('M')
    loc2 = loc2.add_suffix('W')
    df = pd.concat([loc1, loc2, loc3], axis=1)
    return(df, loc1, loc2, loc3)

df3, loc1, loc2, loc3 = missingArea(df3)


def ratioAreas(df, loc1, loc2, loc3):
    loc2.columns = loc2.columns.str.replace('W', '')
    loc3.columns = loc3.columns.str.replace('M', '')
    loc4 = loc3/loc1
    loc4 = np.log10(loc4)
    loc4 = loc4.add_suffix('R_log')
    loc4 = loc4.replace(np.inf, np.NaN)
    loc4 = loc4.replace(np.NINF, np.NaN)
    loc4['zeros'] = loc4.isna().sum(axis=1)
    df = pd.concat([df, loc4], axis=1)
    return(df, loc4)

df3, loc4 = ratioAreas(df3, loc1, loc2, loc3)

print(df3)
print(df3.columns)

def qaqc1(df):
    df['W_average'] = np.mean(df.loc[:, 'pos-octlW':'pos_undecW'], axis=1)
    df['W_std'] = np.std(df.loc[:, 'pos-octlW':'pos_undecW'], axis=1)
    df['W_Rstd'] = df['W_std']/df['W_average']*100
    df['W_Rstd_flag'] = np.where(df['W_Rstd'] <= 15, 1, 0)
    
    return(df)

df3 = qaqc1(df3)

def mergEntact(df1, df2, df3, entact):
    # Merge datasets to focus only on masses that are known to be present
    df1['mass_round'] = np.round(df1['Average Mz'], 0)
    entact['mass_round'] = np.round(entact['Monoisotopic_Mass'], 0)
    df = pd.concat([df1, df2, df3, df4], axis=1) 
    dfe = pd.merge(df, entact, on='mass_round')
    dfe['Mass_diff'] = ((np.absolute(dfe['Average Mz'] - dfe['Monoisotopic_Mass']))/dfe['Monoisotopic_Mass'])*10**6
    dfe['Mass_flag'] = np.where(dfe['Mass_diff'] > 15, 1, 0)
    dfe = dfe[dfe['Mass_flag'] == 0]
    return(dfe)

dfe = mergEntact(df1, df2, df3, entact)

print(dfe)
x = dfe['Mass_diff']
sns.distplot(x)

dfe.to_csv('1stExpMergedCleanpos1.0_0.1min_0.00145Da.csv')
























































'''
# Correct areas for areas found in the blank samples
df2['BlankCorr'] = df2.loc[:, 'pos_no injection_1':'pos_no injection_4'].mean(axis=1)

df3 = df3.sub(df2['BlankCorr'], axis=0)
df2 = df2.sub(df2['BlankCorr'], axis=0)

df3[df3 < 0] = 0
df2[df2 < 0] = 0

df3['zeros'] = (df3 == 0).astype(int).sum(axis=1)
'''










