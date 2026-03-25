# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 11:52:37 2023

@author: scamaraditpinto
"""
import Virtual_Lobule_Script_Ansys_V6 as sim
import Read_Results_V4 as Read
import pandas as pd
import os


## Set the fluid simulation pressure input and output
inletPressure = '800 [Pa]'
outletPressure = '500 [Pa]'
## Set the number of second contained in every iteration (!!!!! needs to change the fluent case file if changed !!!!!) 
timeIteration = 1200

## Open input file that contains all the patients you want to run
csvfile = pd.read_excel('InputValues.xlsx').values
# csvfile = pd.read_excel('InputValues_Temp.xlsx').values

## For each patient run a simulation and generate the results of this simulation
for c in range(2, len(csvfile[0])):
    patientNum = str(int(csvfile[0][c]))
    totalIterations = int(csvfile[1][c] * 72)
    totalMass = csvfile[23][c]
    concentration = ['0.0' for x in range(totalIterations)]
    i=1
    
    if totalIterations > 21:
        for r in range (2, 23):
            concentration[i]=str(csvfile[r][c])
            i=i+1
    else:
        for r in range (2, totalIterations+1):
            concentration[i]=str(csvfile[r][c])
            i=i+1

## Set the result directory and create it automatically if not already existing
    directory = 'Results-' + patientNum 
    os.mkdir(directory) 
    
## Run the simulation
    sim.run(0, totalIterations, concentration, inletPressure, outletPressure, timeIteration, directory, False)
## Generate the results
    Read.Print_Results(directory, totalIterations)