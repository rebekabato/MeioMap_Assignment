#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 20:56:05 2019

@author: bebi
"""
# import moduls and functions
import os
import pandas as pd
import PhasingFunction as ph

# set working directory: Please insert the pathway of your working directory, where you keep the files of this program, between the quotes!
os.chdir("")

# load data
readInData = pd.read_csv("input_data_trios.txt", sep='\\t', engine='python')
clearedData = ph.clearData(readInData)

# set egg in first trio as reference
refColName = '1.egg.GType'

# phase cells by reference, save BED files on computer
cellNames = ['1.PB1.GType', '1.PB2.GType', '2.PB1.GType', '2.PB2.GType', '2.egg.GType', '3.PB1.GType', '3.PB2.GType', '3.egg.GType', '4.PB1.GType', '4.PB2.GType', '4.egg.GType']
for cellname in cellNames:
    phases = ph.phase(clearedData, cellname, refColName)
    subset = ph.createSubset(clearedData, phases, refColName)
    aggreagte = ph.aggregateRows(subset)
    save = ph.saveFile(aggreagte, cellname)

