#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 10:28:12 2022

@author: dabrahamsson
"""

import pandas as pd
import numpy as np

df = pd.read_csv('area_neg.csv')

# clean column names

loc1 =  df.loc[:,'n-octane_injectionA':'dichloromethane-W_injectionB']

loc1.columns = loc1.columns.str.replace('_injectionA', '_1')
loc1.columns = loc1.columns.str.replace('_injectionB', '_2')

loc1.columns = loc1.columns.str.replace('1-octanol'      , '1-octl')
loc1.columns = loc1.columns.str.replace('butyl_acetate'  , 'butac')
loc1.columns = loc1.columns.str.replace('chloroform'     , 'chlor')
loc1.columns = loc1.columns.str.replace('cyclohexane'    , 'chx')
loc1.columns = loc1.columns.str.replace('dichloromethane', 'dcm')
loc1.columns = loc1.columns.str.replace('triolein'       , 'glyctriol')
loc1.columns = loc1.columns.str.replace('n-hexane'       , 'hx')
loc1.columns = loc1.columns.str.replace('n-octane'       , 'octn')
loc1.columns = loc1.columns.str.replace('oleyl_alcohol'  , 'olalc')
loc1.columns = loc1.columns.str.replace('toluene'        , 'tol')
loc1.columns = loc1.columns.str.replace('n-undecane'     , 'undec')

loc1.columns = loc1.columns.str.replace('-W', 'W')
loc1.columns = 'neg_' + loc1.columns

loc1 = loc1[['neg_1-octl_1','neg_1-octl_2','neg_1-octlW_1','neg_1-octlW_2',	
             'neg_butac_1', 'neg_butac_2',	'neg_butacW_1', 'neg_butacW_2', 	
             'neg_chlor_1', 'neg_chlor_2',	'neg_chlorW_1','neg_chlorW_2',	
             'neg_chx_1',	'neg_chx_2',	'neg_chxW_1',	'neg_chxW_2',	'neg_dcm_1',	
             'neg_dcm_2',	'neg_dcmW_1',	'neg_dcmW_2',	'neg_glyctriol_1',	
             'neg_glyctriol_2', 'neg_glyctriolW_1','neg_glyctriolW_2',	
             'neg_hx_1', 'neg_hx_2',	'neg_hxW_1',	'neg_hxW_2',	'neg_octn_1',	
             'neg_octn_2',	'neg_octnW_1',	'neg_octnW_2',	'neg_olalc_1',	'neg_olalc_2',	
             'neg_olalcW_1',	'neg_olalcW_2',	'neg_tol_1',	'neg_tol_2',	
             'neg_tolW_1',	'neg_tolW_2',	'neg_undec_1',	'neg_undec_2',	'neg_undecW_1',	'neg_undecW_2']]

print(loc1.columns)

df =  pd.concat([df, loc1], axis=1)

df.to_csv('area_neg_ready.csv')
