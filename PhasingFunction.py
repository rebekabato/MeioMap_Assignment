#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  2 23:15:09 2019

@author: bebi
"""
import pandas as pd
import numpy as np


# clearData() useage: clean raw data - eliminate non-significant elements
# input: loaded data
# output: data frame containing only significant elements
def clearData(data):
    # keep rows containing informative heterozygous maternal SNP
    data = data[data['gDNA.GType'] == 'AB']
    # eliminate rows containing NCs and heterozygous genotype in reference column
    data = data[data['1.egg.GType'] != 'NC']
    data = data[data['1.egg.GType'] != 'AB']
    return data



# pase() useage: phase cells by reference
# input: data frame, name of the column of the haplotype being examined, reference
# output: np.array of the result after phasing
def phase(clearedData, colName, refColName):
    result = np.array([])
    for index, row in clearedData.iterrows():
        if row[colName] == row[refColName]:
            result = np.append(result, 1)
            
        elif row[colName] == "AB":
            result = np.append(result, 0.5)
            
        elif row[colName] == "NC":
            result = np.append(result, "NC")
        
        else:
            result = np.append(result, 0)
            
    return result



# createSubset() useage: add the array of the result of the phasing process as a column to the new subset of the data frame
# input: data frame, phased column of the cell being examined, name of the column of the haplotype being examined
# output: subset for a particular cell - Chromosome, Position, Phase mark
def createSubset(clearedData, phasedColumn, colName):
    subset = clearedData[['Chr', 'Position']]
    subset['Phase'] = phasedColumn
    return subset



# aggregateRows() useage: create aggregated data frame for a particular subset,
#   aggregation of rows by Phase column, extract the starting and the ending 
#   position of the chromosome
# input: subset of a particular cell - Chromosome, Position, Phase mark
# output: data frame with aggregated rows by Phase marks
def aggregateRows(subset):
    startingPositions = np.array([])
    endPositions = np.array([])
    phases = np.array([])
    lastPhase = None
    currentEndPosition = None
    for index, row in subset.iterrows():
        if lastPhase == None:
            currentEndPosition = row['Position']
            lastPhase = row['Phase']
            phases = np.append(phases, lastPhase)
            startingPositions = np.append(startingPositions, str(row['Position']))
        else:
            if lastPhase != row['Phase']:
                lastPhase = row['Phase']
                startingPositions = np.append(startingPositions, row['Position'])
                phases = np.append(phases, lastPhase)
                endPositions = np.append(endPositions, str(currentEndPosition))
        currentEndPosition = row['Position']
    endPositions = np.append(endPositions, str(currentEndPosition))
    ch1Col = np.full(np.size(startingPositions), 'chr1')
    matrix = np.column_stack((ch1Col, startingPositions, endPositions, phases))
    aggregate = pd.DataFrame(matrix, columns=['chr', 'start', 'end', 'phase'])
    return aggregate



# saveFile() useage: save aggregated dataframe to tab-delimited BED file
# input: aggregated data frame of a particular cell
# output: saved BED file on computer
def saveFile(aggregate, colName):
    aggregate.to_csv(colName + '.bed', sep="\t")
    

