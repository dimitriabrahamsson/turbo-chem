#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 18:19:13 2020

@author: dabrahamsson
"""

#from __future__ import absolute_import, division, print_function

import tensorflow as tf
from tensorflow import keras
#import sklearn
#from sklearn.metrics import mean_absolute_error
print(tf.__version__)

# Import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import RobustScaler

df = pd.read_csv('TurboChemDB5.2.csv')
loc1 = df.loc[:, 'Query':'MOLECULAR_FORMULA']
loc2 = df.loc[:, 'Al_COO':'urea']
loc2 = loc2[['Al_COO', 'Ar_COO', 'Ar_N', 'Ar_NH', 'Ar_OH', 'COO', 'COO2', 'C_O', 
             'C_O_noCOO', 'NH0', 'NH1','NH2', 'Nhpyrrole', 'alkyl_halide', 'amide', 'benzene',
             'bicyclic', 'ester', 'ether', 'halogen', 'phenol', 'phenol_noOrthoHbond']]
loc3 = df.loc[:, 'H':'Br']
loc4 = df.loc[:, 'i_Al_COO':'i_urea']
df = pd.concat([loc1, loc2, loc3, loc4], axis = 1)

def read_and_clean():
    # Read datasets
    dataset = df
    dataset = dataset.loc[:, ~dataset.columns.str.contains('%')] #remove methanol mixtures
    # Show dataframe column names
    dataset.columns.values
    # Remove inorganics and charged molecules
    dataset = dataset[dataset['QSAR_READY_SMILES'].str.contains('CC')]
    dataset = dataset[~dataset['QSAR_READY_SMILES'].str.contains('\+')]
    dataset = dataset[~dataset['QSAR_READY_SMILES'].str.contains('\-')]
    return (dataset)

dataset = read_and_clean()
dataset.columns.values

dataset.to_csv('clean_dataset.csv')

def prep_X(dataset):
    # Assign X and y variables
    X1 = dataset[['1-octanol', 'butyl acetate', 'chloroform', 'cyclohexane',
                  'dichloromethane', 'triolein', 'n-hexane', 'n-octane', 'oleylalcohol',
                  'toluene', 'n-undecane']] # Extract solvent systems from dataframe
    
    sol = ['1-octanol', 'butyl acetate', 'chloroform', 'cyclohexane',
           'dichloromethane', 'triolein', 'n-hexane', 'n-octane', 'oleylalcohol',
           'toluene', 'n-undecane']
    X1.to_csv('sample_solvents_i1.csv') # Print out selection   
    
    scaler = StandardScaler()
    X1 = X1.T
    scaler.fit(X1)
    X1 = scaler.transform(X1)
    X1 = X1.T
    X1 = pd.DataFrame(X1, columns = sol)
    return(X1)

X1 = prep_X(dataset)
    
def assign_X(X1, dataset):    
    X1 = X1.reset_index(drop=True)
    X2 = dataset.loc[:, 'H':'Br'] # Select atoms
    X2 = X2.reset_index(drop=True)
    X3 = dataset.loc[:, 'MONOISOTOPIC_MASS'] # Select molecular mass
    X3 = X3.reset_index(drop=True)
    X4 = dataset.loc[:, 'INCHIKEY'] #  Select InChiKey
    X4 = X4.reset_index(drop=True)
    X5 = dataset.loc[:, 'i_Al_COO':'i_urea']
    X5 = X5.reset_index(drop=True)
    X = pd.concat([X1, X2, X3, X4, X5], axis=1) # Compose a dataframe that contains the X parameters
    return (X)

X = assign_X(X1, dataset)
print(X.columns)
    
def assign_y():
    y1 = dataset.loc[:, 'Al_COO':'phenol_noOrthoHbond'] # Extract RDKit fragments from dataframe
    y2 = dataset.loc[:, 'INCHIKEY']
    y = pd.concat([y1, y2], axis=1) # Compose a dataframe that contains the y parameters
    return (y)
    
y = assign_y()


def split():
    # Split dataframe into training and testing set
    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=None)
    
    y_tr_in = y_train.loc[:, 'INCHIKEY'].reset_index() # Save indexes for printing results later
    y_ts_in = y_test.loc[:, 'INCHIKEY'].reset_index()
    
    X_train = X_train.set_index('INCHIKEY')
    X_test = X_test.set_index('INCHIKEY')
    y_train = y_train.set_index('INCHIKEY')
    y_test= y_test.set_index('INCHIKEY')
    
    X_train_0 = X_train.values
    X_test_0 = X_test.values
    y_train = y_train.values
    y_test = y_test.values
    
    return (X_train_0, X_test_0, y_train, y_test, y_tr_in, y_ts_in)

X_train_0, X_test_0, y_train, y_test, y_tr_in, y_ts_in = split()



# Read experimental data
exp = pd.read_excel('ExpDataFragment_atoms2.0.xlsx')
exp.columns = exp.columns.str.replace('Indigo_Inchi_Key', 'INCHIKEY')
exp.columns = exp.columns.str.replace('Monoisotopic_Mass', 'MONOISOTOPIC_MASS')
exp.columns = exp.columns.str.replace('_exp_std', '')
exp = exp.loc[:, exp.columns.isin(X.columns)]
exp = exp.reindex(X.columns, axis=1)
exp_in = exp.loc[:, 'INCHIKEY'].reset_index()
exp = exp.set_index('INCHIKEY')
print(exp.columns)
exp.to_csv('ExpDataFragAtomsR.csv')
exp = exp.values


# Normalize data
def normalize():
    sc = StandardScaler()
    X_train = sc.fit_transform(X_train_0)
    X_test = sc.transform(X_test_0)
    exp_test = sc.transform(exp)
    return(X_train, X_test, exp_test)

X_train, X_test, exp_test = normalize()

# Compile the ANN model
def build_model():
  model = keras.Sequential([
    keras.layers.Flatten(input_shape=(X_train.shape[1],)),
    keras.layers.Dense(500, activation='relu'),
    keras.layers.Dense(500, activation='relu'),
    keras.layers.Dense(500, activation='relu'),
    keras.layers.Dense(500, activation='relu'),
    keras.layers.Dense(500, activation='relu'),
    keras.layers.Dense(500, activation='relu'),
    keras.layers.Dense(500, activation='relu'),
    keras.layers.Dense(500, activation='relu'),
    keras.layers.Dense(500, activation='relu'),
    keras.layers.Dense(500, activation='relu'),
    keras.layers.Dropout(0.2),
    keras.layers.Dense(500, activation='exponential'),
    keras.layers.Dense(22)
  ])

  optimizer = tf.optimizers.Adamax(0.001)

  model.compile(loss='mae',
                optimizer=optimizer,
                metrics=['mae'])
  return model

model = build_model()
model.summary()


# Display training progress by printing a single dot for each completed epoch
class PrintDot(keras.callbacks.Callback):
  def on_epoch_end(self, epoch, logs):
    if epoch % 100 == 0: print('')
    print('.', end='')

EPOCHS = 100

# Store training stats
history = model.fit(X_train, y_train, epochs=EPOCHS,
                    validation_split=0.2, verbose=0,
                    callbacks=[PrintDot()])

# Show a plot with the errors for the training and testing set
fig = plt.figure()
ax = plt.subplot()
ax.plot(history.epoch, np.array(history.history['mae']),
           label='Train Loss')
ax.plot(history.epoch, np.array(history.history['val_mae']),
           label = 'Val Loss')
plt.xlabel('Epoch')
plt.ylabel('Mean Abs Error')
plt.legend()
plt.ylim([0, 2.5])
plt.xlim([0, 100])
plt.show()
#fig.savefig('train_val_entact_500_20.png', dpi=400)

def print_mae():
    [loss, mae] = model.evaluate(X_test, y_test, verbose=0)
    print("Testing set Mean Abs Error: {:7.2f}".format(mae * 1))
    [loss, mae] = model.evaluate(X_train, y_train, verbose=0)
    print("Training set Mean Abs Error: {:7.2f}".format(mae * 1))
    return()

print_mae()

# Use model to predict data for training set
def predict_train():
    y_pred = model.predict(X_train)
    y_names = dataset.loc[:, 'Al_COO':'phenol_noOrthoHbond'].columns.values
    df_y_train = pd.DataFrame(y_train, columns=y_names)
    df_y_pred1 = pd.DataFrame(y_pred, columns=y_names)
    df_y_train = df_y_train.add_suffix('_x')
    df_y_pred1 = df_y_pred1.add_suffix('_y')
    df_y1 = pd.concat([df_y_train, df_y_pred1], axis=1)
    return(df_y1, y_pred, y_names)

df_y1, y_pred, y_names = predict_train()


def graphs(df_y1, y_train, y_pred):
    import seaborn as sns
    x = df_y1['Al_COO_x']
    y = df_y1['Al_COO_y']
    ax = sns.scatterplot(x=x, y=y, data=df_y1, color='dodgerblue', alpha=0.1, s=80)
    _ = plt.plot([0, 25], [0, 25])
    plt.show()
    
    x = df_y1['benzene_x']
    y = df_y1['benzene_y']
    ax = sns.scatterplot(x=x, y=y, data=df_y1, color='dodgerblue', alpha=0.1, s=80)
    _ = plt.plot([0, 10], [0, 10])
    plt.show()
    
    x = df_y1['ether_x']
    y = df_y1['ether_y']
    ax = sns.scatterplot(x=x, y=y, data=df_y1, color='dodgerblue', alpha=0.1, s=80)
    _ = plt.plot([0, 10], [0, 10])
    plt.show()
    
    error = y_pred - y_train
    plt.hist(error, bins = 10)
    plt.xlabel("Prediction Error")
    _ = plt.ylabel("Count")
    plt.show()
    return()

graphs(df_y1, y_train, y_pred)
    
# Use model to make predictions for the testing set
    
def predict_test():
    y_pred = model.predict(X_test)
    df_y_test = pd.DataFrame(y_test, columns=y_names)
    df_y_pred2 = pd.DataFrame(y_pred, columns=y_names)
    df_y_test = df_y_test.add_suffix('_x')
    df_y_pred2 = df_y_pred2.add_suffix('_y')
    df_y2 = pd.concat([df_y_test, df_y_pred2], axis=1)
    return(y_pred, df_y2)

y_pred, df_y2 = predict_test()

graphs(df_y2, y_test, y_pred)


# Make predictions for experimental data
exp_pred = model.predict(exp_test)
df_exp = pd.DataFrame(exp_pred, columns=y_names)
df_exp.to_csv('PredsExpData5.csv')


def print_results(df_y1, df_y2):
    df_y1 = pd.concat([df_y1, y_tr_in], axis=1)
    df_y2 = pd.concat([df_y2, y_ts_in], axis=1)
    df_y1 = df_y1.set_index('INCHIKEY')
    df_y2 = df_y2.set_index('INCHIKEY')
    df_y1.to_csv('PredsTrain5.csv')
    df_y2.to_csv('PredsTest5.csv')
    return ()

print_results(df_y1, df_y2)








