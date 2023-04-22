#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 13:24:25 2022

@author: dabrahamsson
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import linear_model
from random import sample

for i in range(0,100):

    def read_and_clean():
        df = pd.read_csv('BloodExposomeLSERfragment_atoms_R.csv')
        df = df[['wet 1-octanol', 'wet butyl acetate', 'CHCl3', 
                  'wet & dry Cyclohexane', 'wet & dry Dichloromethane',
                 'wet triolein', 'n-hexane', 'n-octane', 
                 'wet & dry Oleylalcohol', 'wet & dry Toluene', 'n-undecane']]
        
        df.columns = df.columns.str.replace('wet & dry ', '')
        df.columns = df.columns.str.replace('wet ', '')
        df.columns = df.columns.str.replace('CHCl3', 'chloroform')
        df.columns = df.columns.str.lower()
        sol = df.columns
        print(df.columns)
        
        dfR = pd.read_csv('1stExpMergedCleanNeg.csv')
        dfR = dfR.replace(np.inf, np.NaN)
        dfR['chem_id'] = np.round(dfR['Average Mz'], 3).astype(str) + '@' + np.round(dfR['Average Rt(min)'], 3).astype(str)
        dfR = dfR.set_index('chem_id')
        dfR = dfR.loc[:, 'neg-octlR_log':'neg_undecR_log']
        dfR.columns = sol
        dfR = dfR.reset_index()
        print(dfR.columns)
        return(df, dfR, sol)
    
    df, dfR, sol = read_and_clean()
    
    solS1 = sample(sorted(sol), 3)
    solS2 = sol[~sol.isin(solS1)]
    
    s = ''
    for i in solS1:
        s+=str(i[:4])
        
    print(s)
    print(solS1)
    print(solS2)
    
    
    def read_data(df, dfR, solS1, solS2):
        dfR1 = dfR.dropna(subset=solS1)
        X = df[solS1].values
        y = df[solS2].values
        X1 = dfR1[solS1].values 
        chemid = pd.DataFrame(dfR1['chem_id'], columns=['chem_id'])
        chemid = chemid.reset_index(drop=True)
        return(X, y, X1, chemid)
    
    
    def predict(X, y, X1, s4):
        model_ols =  linear_model.LinearRegression(normalize=True)
        model_ols.fit(X,y)
        pred = pd.DataFrame(model_ols.predict(X), columns=['predicted_{}'.format(s4)]) # Create new dataframe of column'Predicted Price'
        actual = pd.DataFrame(y, columns=['true_{}'.format(s4)])
        actual = actual.reset_index(drop=True) # Drop the index so that we can concat it, to create new dataframe
        comp = pd.concat([actual,pred],axis =1)
        plt.scatter(y, model_ols.predict(X))
        
        # Make predictions for experimental data
        y1 = model_ols.predict(X1)
        pred = pd.DataFrame(y1, columns=['predicted_{}'.format(s4)]) # Create new dataframe of column'Predicted Price'
        comp1 = pd.concat([chemid,pred], axis=1)
        return(comp1)
    
    for i in solS2:
        X, y, X1, chemid = read_data(df, dfR, solS1, i)
        comp1 = predict(X, y, X1, i)
        comp1.to_csv('{}_'.format(i) + s + '.csv')
        print(comp1)





