#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 18:19:13 2020

@author: dabrahamsson
"""

from __future__ import absolute_import, division, print_function

import tensorflow as tf
from tensorflow import keras
import sklearn
from sklearn.metrics import mean_absolute_error
print(tf.__version__)

# Import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Read datasets
dataset = pd.read_csv('BloodExposomeLSERfragment_atoms_R.csv')
#dataset = dataset.select_dtypes([np.number])

dataset = dataset.loc[:, ~dataset.columns.str.contains('%')] #remove methanol mixtures

# Show dataframe column names
dataset.columns.values

#dataset = dataset[dataset['L_LSER'] == dataset['L_LSER']]
dataset = dataset[dataset['QSAR_READY_SMILES'].str.contains('CC')]
dataset = dataset[~dataset['QSAR_READY_SMILES'].str.contains('\+')]
dataset = dataset[~dataset['QSAR_READY_SMILES'].str.contains('\-')]

# Assign X and y variables
X1 = dataset.loc[:, 'dry triolein':'Octan-1-ol/propylenecarbonate']
X1 = X1.sample(n=10, axis='columns')
X1.to_csv('sample_solvents_i5.csv')
#X1 = pd.read_csv('sample_solvents_i5.csv')
#X1 = X1.drop('Unnamed: 0', axis=1)
X1 = X1.reset_index(drop=True)

X2 = dataset.loc[:, 'H':'Br']
X2 = X2.reset_index(drop=True)

X3 = dataset.loc[:, 'AVERAGE_MASS']
X3 = X3.reset_index(drop=True)

X4 = dataset.loc[:, 'INCHIKEY']
X4 = X4.reset_index(drop=True)

X = pd.concat([X1, X2, X3, X4], axis=1)

y1 = dataset.loc[:, 'Al_COO':'phenol_noOrthoHbond']
y2 = dataset.loc[:, 'INCHIKEY']
y = pd.concat([y1, y2], axis=1)

# Split dataframe into training and testing set
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=None)

X_tr_in = X_train.loc[:, 'INCHIKEY'].reset_index()
X_ts_in = X_test.loc[:, 'INCHIKEY'].reset_index()
y_tr_in = y_train.loc[:, 'INCHIKEY'].reset_index()
y_ts_in = y_test.loc[:, 'INCHIKEY'].reset_index()

X_train = X_train.set_index('INCHIKEY')
X_test = X_test.set_index('INCHIKEY')
y_train = y_train.set_index('INCHIKEY')
y_test= y_test.set_index('INCHIKEY')

X_train = X_train.values
X_test = X_test.values
y_train = y_train.values
y_test = y_test.values

# Normalize data
from sklearn.preprocessing import StandardScaler

sc = StandardScaler()
X_train = sc.fit_transform(X_train)
X_test = sc.transform(X_test)

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

EPOCHS = 200

# Store training stats
history = model.fit(X_train, y_train, epochs=EPOCHS,
                    validation_split=0.2, verbose=0,
                    callbacks=[PrintDot()])


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
plt.xlim([0, 200])
plt.show()
#fig.savefig('train_val_entact_500_20.png', dpi=400)


'''
# The patience parameter is the amount of epochs to check for improvement
early_stop = keras.callbacks.EarlyStopping(monitor='val_loss', patience=500)

history = model.fit(X_train, y_train, epochs=EPOCHS,
                    validation_split=0.2, verbose=0,
                    callbacks=[early_stop, PrintDot()])

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
plt.show()
#fig.savefig('train_val_entact_20.png', dpi=400)
'''

[loss, mae] = model.evaluate(X_test, y_test, verbose=0)
print("Testing set Mean Abs Error: {:7.2f}".format(mae * 1))

[loss, mae] = model.evaluate(X_train, y_train, verbose=0)
print("Training set Mean Abs Error: {:7.2f}".format(mae * 1))


# Use model to predict data for training set
y_pred = model.predict(X_train)

y_names = dataset.loc[:, 'Al_COO':'phenol_noOrthoHbond'].columns.values

df_y_train = pd.DataFrame(y_train, columns=y_names)
df_y_pred1 = pd.DataFrame(y_pred, columns=y_names)

df_y_train = df_y_train.add_suffix('_x')
df_y_pred1 = df_y_pred1.add_suffix('_y')

df_y1 = pd.concat([df_y_train, df_y_pred1], axis=1)

import seaborn as sns

x = df_y1['Al_OH_x']
y = df_y1['Al_OH_y']
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

# Use model to make predictions for the testing set
y_pred = model.predict(X_test)


df_y_test = pd.DataFrame(y_test, columns=y_names)
df_y_pred2 = pd.DataFrame(y_pred, columns=y_names)


df_y_test = df_y_test.add_suffix('_x')
df_y_pred2 = df_y_pred2.add_suffix('_y')

df_y2 = pd.concat([df_y_test, df_y_pred2], axis=1)

x = df_y2['Al_OH_x']
y = df_y2['Al_OH_y']
ax = sns.scatterplot(x=x, y=y, data=df_y2, color='dodgerblue', alpha=0.1, s=80)
_ = plt.plot([0, 10], [0, 10])
plt.show()

x = df_y2['benzene_x']
y = df_y2['benzene_y']
ax = sns.scatterplot(x=x, y=y, data=df_y2, color='dodgerblue', alpha=0.1, s=80)
_ = plt.plot([0, 10], [0, 10])
plt.show()

x = df_y2['ether_x']
y = df_y2['ether_y']
ax = sns.scatterplot(x=x, y=y, data=df_y2, color='dodgerblue', alpha=0.1, s=80)
_ = plt.plot([0, 10], [0, 10])
plt.show()

error = y_pred - y_test
plt.hist(error, bins = 10)
plt.xlabel("Prediction Error")
_ = plt.ylabel("Count")
plt.show()

df_y1 = pd.concat([df_y1, y_tr_in], axis=1)
df_y2 = pd.concat([df_y2, y_ts_in], axis=1)

df_y1 = df_y1.set_index('INCHIKEY')
df_y2 = df_y2.set_index('INCHIKEY')

df_y1.to_csv('training_results_exposome_i5.4.csv')
df_y2.to_csv('test_results_exposome_i5.4.csv')








