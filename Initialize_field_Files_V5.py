# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 11:26:19 2023

@author: scamaraditpinto
"""
import numpy as np
import shutil
import random as rd


def Initialize_Hepatocytes_File():
    HepatocytesField = open('Field_Hepatocytes.txt' , 'w')
    GSHField = open('Field_Zonation.txt' , 'r+')
    GSHNodes = GSHField.readlines()
    
    for line in GSHNodes:
        GSHMass = line.split(',')
    
        HepatocytesField.write( str(GSHMass[0].rstrip()) + ',' 
                               + '0.0' + ',' 
                               + '0.0' + ',' 
                               + str(GSHMass[1].rstrip()) + ',' 
                               + '1' + ','
                               + '0' + ',' 
                               + '0' + ',' 
                               + '0' + ',' 
                               + '12' + ',' 
                               + '9' + ','
                               + '0' + ','
                               + str(GSHMass[2].rstrip()) + ','
                               + str(GSHMass[3].rstrip()) + ','
                               + str(GSHMass[4].rstrip()) + ','
                               + str(GSHMass[5].rstrip()) + '\n')
        
    HepatocytesField.close()
    GSHField.close()
    return

def Initialize_Zonation_File():
    zonationField = open('Field_Zonation.txt', 'w')
    nodeField = open('Field_Nodes.txt', 'r+')
    Nodes = nodeField.readlines()
    
    lobuleR = 0.75 ## mm
    centralveinR = 0.075 ## mm
    ## Less GSH in central          Initial value of the model = 0.7e-14 for depletion at 6g
    valueGSHInitCentral = 0.5e-14 ## mol
    valueGSHInitPortal = 1e-14 ## mol
    ## Less production in central   Initial value of the model = 1.375e-14
    # percentVariation = 0.2
    valueGSHProcducRateCentral = 1.175e-14/8.64e4# - (1.175e-14/8.64e4 * percentVariation) ## mol
    valueGSHProductRatePortal = 1.575e-14/8.64e4# + (1.175e-14/8.64e4 * percentVariation) ## mol
    ## Less binding in portal       Initial value of the model = 1.6e18
    # percentVariation = 0.2
    valueGSHBindingRateCentral = 1.4e18/8.64e4#  - (1.6e18/8.64e4 * percentVariation) ## mol
    valueGSHBindingRatePortal = 1.8e18/8.64e4# + (1.6e18/8.64e4 * percentVariation) ## mol
    ## More Regen in Portal         Initial value of the model = 1
    # percentVariation = 0.2
    valueRegenRateCentral = 0.9/8.64e4# - (1/8.64e4 * percentVariation) ## mol
    valueRegenRatePortal = 1.1/8.64e4# + (1/8.64e4 * percentVariation)## mol 
    
    a = (valueGSHInitPortal-valueGSHInitCentral) / ((lobuleR*1e-3)-(centralveinR*1e-3))
    b = valueGSHInitCentral - (a*centralveinR*1e-3)
    c = (valueGSHProductRatePortal-valueGSHProcducRateCentral) / ((lobuleR*1e-3)-(centralveinR*1e-3))
    d = valueGSHProcducRateCentral - (c*centralveinR*1e-3)
    e = (valueGSHBindingRatePortal-valueGSHBindingRateCentral) / ((lobuleR*1e-3)-(centralveinR*1e-3))
    f = valueGSHBindingRateCentral - (e*centralveinR*1e-3)
    g = (valueRegenRatePortal-valueRegenRateCentral) / ((lobuleR*1e-3)-(centralveinR*1e-3))
    h = valueRegenRateCentral - (g*centralveinR*1e-3)
    
    for line in Nodes:
        nodeCoords = line.split(',')
        GSHInitValue = a * np.sqrt((float(nodeCoords[1].rstrip()) * float(nodeCoords[1].rstrip())) + (float(nodeCoords[2].rstrip()) * float(nodeCoords[2].rstrip()))) + b
        GSHProducRateValue = c * np.sqrt((float(nodeCoords[1].rstrip()) * float(nodeCoords[1].rstrip())) + (float(nodeCoords[2].rstrip()) * float(nodeCoords[2].rstrip()))) + d
        GSHBindingRateValue = e * np.sqrt((float(nodeCoords[1].rstrip()) * float(nodeCoords[1].rstrip())) + (float(nodeCoords[2].rstrip()) * float(nodeCoords[2].rstrip()))) + f
        regenRateValue = g * np.sqrt((float(nodeCoords[1].rstrip()) * float(nodeCoords[1].rstrip())) + (float(nodeCoords[2].rstrip()) * float(nodeCoords[2].rstrip()))) + h
        random = rd.random()
        
        zonationField.write( nodeCoords[0] + ',' + str(GSHInitValue) + ',' + str(GSHProducRateValue) + ',' + str(GSHBindingRateValue) + ',' + str(regenRateValue) + ',' + str(random) + '\n')
        
    zonationField.close()
    nodeField.close()
    return

def Initialize_Porosity_File():
    cell_spericity = 0.95
    cell_diameter = 25e-6
    porosity = 0.143
    
    inertial_res = 2*1.75*(1-porosity) / (cell_spericity * cell_diameter * np.power(porosity,3))
    # print (inertial_res)
    viscous_res = 150*np.power(1-porosity,2) / (np.power(cell_spericity,2) * np.power(cell_diameter,2) * np.power(porosity,3))
    # print (viscous_res)
    shutil.copyfile("Profile_porosity-Template.prof", 'Profile_porosity.prof')
    porosityField = open('Profile_porosity.prof' , 'a+')
    nodeField = open('Field_Nodes.txt' , 'r+')
    Nodes = nodeField.readlines()
    
    porosityField.write('(porosity \n')
    for line in Nodes:
        porosityField.write(str(0.143) + '\n')
    porosityField.write(') \n')
    
    porosityField.write('(Inertial-resistance \n')

    for line in Nodes:
        porosityField.write(str(inertial_res) + '\n')
    porosityField.write(') \n')
    
    porosityField.write('(viscous-resistance \n')
    
    for line in Nodes:
        porosityField.write(str(viscous_res) + '\n')
    porosityField.write(')) \n')
    
    porosityField.close()
    nodeField.close()
    return

def Initialize_Status_File():
    statusFile = open('hepatocytes_status_Hepatocytes.txt', 'w')
    nodeField = open('Field_Nodes.txt', 'r+')
    Nodes = nodeField.readlines()
    for line in Nodes:
        statusFile.write('0.0\n')
    statusFile.close()
    nodeField.close()
    return