#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 10:28:12 2022

@author: dabrahamsson
"""

import pandas as pd
import numpy as np

df = pd.read_csv('area_pos.csv')

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
loc1.columns = 'pos_' + loc1.columns

loc1 = loc1[['pos_1-octl_1','pos_1-octl_2','pos_1-octlW_1','pos_1-octlW_2',	
             'pos_butac_1', 'pos_butac_2',	'pos_butacW_1', 'pos_butacW_2', 	
             'pos_chlor_1', 'pos_chlor_2',	'pos_chlorW_1','pos_chlorW_2',	
             'pos_chx_1',	'pos_chx_2',	'pos_chxW_1',	'pos_chxW_2',	'pos_dcm_1',	
             'pos_dcm_2',	'pos_dcmW_1',	'pos_dcmW_2',	'pos_glyctriol_1',	
             'pos_glyctriol_2', 'pos_glyctriolW_1','pos_glyctriolW_2',	
             'pos_hx_1', 'pos_hx_2',	'pos_hxW_1',	'pos_hxW_2',	'pos_octn_1',	
             'pos_octn_2',	'pos_octnW_1',	'pos_octnW_2',	'pos_olalc_1',	'pos_olalc_2',	
             'pos_olalcW_1',	'pos_olalcW_2',	'pos_tol_1',	'pos_tol_2',	
             'pos_tolW_1',	'pos_tolW_2',	'pos_undec_1',	'pos_undec_2',	'pos_undecW_1',	'pos_undecW_2']]

print(loc1.columns)

df =  pd.concat([df, loc1], axis=1)

df.to_csv('area_pos_ready.csv')
