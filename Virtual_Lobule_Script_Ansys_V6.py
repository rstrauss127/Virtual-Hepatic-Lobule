# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 14:46:38 2023

@author: scamaraditpinto

ADDED REGEN, ZONATION, REMOVED NAC and random state
Classification of status is randomized
"""

import ansys.fluent.core as pyfluent
import shutil
import numpy as np
from time import time



# MALD model tools
import MALD_model_Scaled_V5 as MALD
import Initialize_field_Files_V5 as Init
# import Initialize_field_Files_V5_NoZonation as Init


# Profile file path
profile_path = "Profile_porosity-Template.prof"
result_path = 'Lobule-Flow-Drug-Simulation-1200'
case_path = 'Lobule-Flow-Drug-Simulation-v4-1200.cas.h5'

def run (firstIteration, totalIterations, drugConcentration, inletPressure, outletPressure, timeIteration, fileNamePrefixe, runCFD):
    ## Mark time of the start of the simualtion
    time0 = time()
    
    ## porosity values
    porosity_normal = 0.143
    porosity_lysed = 0.246
    cell_spericity = 0.95
    cell_normal_diameter = 25e-6
    cell_damaged_diameter = 25e-6
    
    ## Initialize data files based on node list ##
    Init.Initialize_Zonation_File()
    Init.Initialize_Hepatocytes_File()
    Init.Initialize_Porosity_File()
    Init.Initialize_Status_File()
    
    ## Run loop to iterrate
    for n in range (totalIterations):
## Mark the begining of the iteration
        print('Begining of iteration ' + str(n))
        
        if drugConcentration[n] != '0.0' and runCFD:
## Open ansys solver and case
            solver = pyfluent.launch_fluent(precision="double", version="2d", processor_count=4, mode="solver", show_gui=False, cleanup_on_exit=True, case_filepath=case_path)
            ## Set up inlet pressure
            solver.setup.named_expressions['parameter_inlet_pressure'].definition = inletPressure
            ## Set up outlet pressure
            solver.setup.named_expressions['parameter_outlet_pressure'].definition = outletPressure
            ## Set up number of total iteration in the flow simulation
            solver.solution.run_calculation.transient_controls.time_step_count = timeIteration
            ## Set up drug percentile of composition of blood
            solver.setup.named_expressions['parameter_inlet_drug'].definition = drugConcentration[n] 
                
## Initiallize and run ##
            solver.solution.initialization.hybrid_initialize()
            solver.solution.run_calculation.calculate()
            solver.exit()
        
## Initialize values for hepatocytes and porosity arrays based on previous values
        ## Open the file that contains the values of the hepatocytes status
        hepatocytesField = open('Field_Hepatocytes.txt' , 'r+')
        ## Separate every line
        hepatocytesLines = hepatocytesField.readlines()
        ## Initialize storage array for value in the hepatocytes
        hepatocytes_array = np.zeros((len(hepatocytesLines), 10))
        hepatocytesZonation_array = np.zeros((len(hepatocytesLines), 3))
        hepatocytesRandom_array = np.zeros(len(hepatocytesLines))
        ## Put every value in every line into the storage array
        for i in range (len(hepatocytesLines)):
            splitLines = hepatocytesLines[i].split(',')
            for x in range (10):
                hepatocytes_array[i][x] = float(splitLines[x+1])
            hepatocytesZonation_array[i][0] = float(splitLines[-4])
            hepatocytesZonation_array[i][1] = float(splitLines[-3])
            hepatocytesZonation_array[i][2] = float(splitLines[-2])
            hepatocytesRandom_array[i] = float(splitLines[-1])
        ## Close the file and save it as used values
        
        concentrationAPAP = np.zeros(len(hepatocytesLines))
        
        hepatocytesField.close()
        shutil.copyfile('Field_Hepatocytes.txt', fileNamePrefixe +'/_Used-Field_Hepatocytes-%d.txt' %(n))
        
        if drugConcentration[n] != '0.0' and runCFD:
## initialize concentration added to each hepatocytes based on the flow analysis results from ansys
            ## Open the result file
            concentrationField = open(result_path , 'r+')
            ## Separate every line 
            concentrationLines = concentrationField.readlines()
            ## Store the results in the storage array 
            for i in range (len(concentrationLines)-1):
                concentrationAPAP[i] = float(concentrationLines[i+1].split(',')[-1])
            ## Close the file and save it as used values
            concentrationField.close()
            shutil.copyfile(result_path, fileNamePrefixe +'/_Used-' + result_path + '-%d.txt' %(n))
            
        elif drugConcentration[n] != '0.0':
            concentrationField = open(fileNamePrefixe +'/_Used-' + result_path + '-%d.txt' %(n) , 'r+')
            concentrationLines = concentrationField.readlines()
            ## Store the results in the storage array 
            for i in range (len(concentrationLines)-1):
                concentrationAPAP[i] = float(concentrationLines[i+1].split(',')[-1])
            ## Close the file and save it as used values
            concentrationField.close()
            
## Add currently added concentration to previous value
        for j in range (len(concentrationAPAP)):
            hepatocytes_array[j][0] = hepatocytes_array[j][0] + (concentrationAPAP[j] *3.4e-15*timeIteration)
    
## initialize porosity values and files to be used in the drug kinetics simulation
        ## Copy profile porosity file to use as template
        shutil.copyfile(profile_path, 'Profile_porosity.prof')
        ## Open the new porosity file created to add values
        porosityField = open('Profile_porosity.prof' , 'a+')        
        ## Initialize storage array to store porosity values of every hepatocytes        
        porosity_array = np.zeros((len(hepatocytesLines), 4))
        ## Open hepatocytes filed file to rewrite values with new results
        hepatocytesField = open('Field_Hepatocytes.txt' , 'w')


## open last state of cells
        hepatocytesStatusFile = open('hepatocytes_status_Hepatocytes.txt' , 'r+')
        statusLine = hepatocytesStatusFile.readlines()
        hepatocytesPreviousStatus_array = np.zeros(len(statusLine))
        for i in range (len(statusLine)):
            hepatocytesPreviousStatus_array[i] = str(statusLine[i].split('\n')[0])
        hepatocytesStatusFile.close()
        
## Processing the drug kinetics simulation for every hepatocytes
        for k in range (len(hepatocytes_array[:])):
            ## Run the model based on the hepatocyte status for the same duration of the flow simulation
            solution = MALD.RunMALD(hepatocytes_array[k], timeIteration, hepatocytesZonation_array[k][0], hepatocytesZonation_array[k][1],hepatocytesZonation_array[k][2])
            ## Recover values of cells status given by the scaled MALD Simulation
            normal=solution.y[3][-1]       
            damaged=solution.y[4][-1]  
            lysed=solution.y[5][-1]  
            regen=solution.y[6][-1]
            
            random = hepatocytesRandom_array[k]
            ## Set hepatocyte porous properties based on cell status values and store into storage array         
            if random < regen :
                ## Set porosity value to normal status as healthy and damaged have the same void ratio
                porosity = porosity_normal
                ## Set the cell diameter for viscousresistance value caculation
                cell_diameter = cell_normal_diameter
                ## Set status for graphical visualization value to healthy = 0
                status = '3.0'
            else:
                total = damaged + lysed + normal
                if random < (lysed/total) and str(hepatocytesPreviousStatus_array[k]) != '3.0':
                    ## Set porosity value to lysed status
                    porosity = porosity_lysed
                    ## Set the cell diameter for viscousresistance value caculation
                    cell_diameter = cell_damaged_diameter
                    ## Set status for graphical visualization value to lysed = 2
                    status = '2.0'
                else:
                    total = normal + damaged
                    if random < (damaged/total) and str(hepatocytesPreviousStatus_array[k]) != '2.0' and str(hepatocytesPreviousStatus_array[k]) != '3.0':
                        ## Set porosity value to normal status as healthy and damaged have the same void ratio
                        porosity = porosity_normal
                        ## Set the cell diameter for viscousresistance value caculation
                        cell_diameter = cell_damaged_diameter
                        ## Set status for graphical visualization value to damaged = 1
                        status = '1.0'
                    else:
                        status = hepatocytesPreviousStatus_array[k]
                        porosity = porosity_normal
                        ## Set the cell diameter for viscousresistance value caculation
                        cell_diameter = cell_damaged_diameter
                        
            ## Set up porosity array value where [0] = void ratio, [1]=viscous resistance, [2]=inertial resistance, [3]=status for graphical visualization
            porosity_array[k][0] = porosity
            porosity_array[k][1] = 2*1.75*(1-porosity) / (cell_spericity * cell_diameter * np.power(porosity,3)) # calculated from ergun equation
            porosity_array[k][2] = 150*np.power(1-porosity,2) / (np.power(cell_spericity,2) * np.power(cell_diameter,2) * np.power(porosity,3)) # calculated from ergun equation
            porosity_array[k][3] = status
            
            ## save new state of the hepatocyte           
            hepatocytesField.write( str(k+1) + ','
                                   + str(solution.y[0][-1]) + ',' 
                                   + str(solution.y[1][-1]) + ',' 
                                   + str(solution.y[2][-1]) + ',' 
                                   + str(solution.y[3][-1]) + ','
                                   + str(solution.y[4][-1]) + ',' 
                                   + str(solution.y[5][-1]) + ',' 
                                   + str(solution.y[6][-1]) + ',' 
                                   + str(solution.y[7][-1]) + ','
                                   + str(solution.y[8][-1]) + ','
                                   + str(solution.y[9][-1]) + ','
                                   + str(hepatocytesZonation_array[k][0]) + ','
                                   + str(hepatocytesZonation_array[k][1]) + ','
                                   + str(hepatocytesZonation_array[k][2]) + ','
                                   + str(hepatocytesRandom_array[k]) + '\n')
## Save results from the drug simulation in porosity field file
        ## Save void ratio value for every hepatocytes
        porosityField.write('(porosity \n')
        for l in range (len(porosity_array[:])):
            porosityField.write(str(porosity_array[l][0]) + '\n')
        porosityField.write(') \n')
        ## Save viscous resistance for every hepatocytes
        porosityField.write('(Inertial-resistance \n')
        for l in range (len(porosity_array[:])):
            porosityField.write(str(porosity_array[l][1]) + '\n')
        porosityField.write(') \n')
        ## Save inertial resistance for every hepatocytes
        porosityField.write('(viscous-resistance \n')
        for l in range (len(porosity_array[:])):
            porosityField.write(str(porosity_array[l][2]) + '\n')
        porosityField.write(')) \n')
## close used files
        porosityField.close()
        ## Copy porosity profile file for saving
        shutil.copyfile('Profile_porosity.prof', fileNamePrefixe +'/_Used-Profile_porosity-%d.prof' %(n))
        hepatocytesField.close()
        
## Create visualization of the drug kinetics results showing the status of every cell
        ## Create and open hepatocyte staus file to strore every hepatocyte status
        hepatocytesStatusFile = open('hepatocytes_status_Hepatocytes.txt' , 'w')
        ## Store status values for every hepatocyte
        for l in range (len(porosity_array[:])):
            hepatocytesStatusFile.write(str(porosity_array[l][3]) + '\n')
        ## Close file
        hepatocytesStatusFile.close()
        ## Copy status file for saving
        shutil.copyfile('hepatocytes_status_Hepatocytes.txt', fileNamePrefixe +'/_hepatocytes_status_Hepatocytes-%d.txt' %(n))
## Mark end of the iteration
        print('End of iteration ' + str(n))
        
## Mark time at the end of the simualtion
    time1 = time()
## Mark end of the simulation and the duration of it
    print ('Total simulation time (s) = ' + str(time1-time0) + 'for simulation : ' + fileNamePrefixe)
    return