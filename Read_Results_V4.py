# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 10:42:08 2023

@author: scamaraditpinto
"""
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import griddata
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

def generate_excel(nodes_file_path, status_files_base_path, output_excel_file_path, num_files):
    wb = Workbook()
    wb.remove(wb.active)  # Erases the default sheet

    # Load nodes data
    nodes = pd.read_csv(nodes_file_path, sep=',', header=None, names=['Hepatocyte', 'X', 'Y'])
    centre = nodes[['X', 'Y']].mean()
    for i in range(1, num_files):
        # Load status data for each iteration
        status_file = f"{status_files_base_path}-{i}.txt"
        status = pd.read_csv(status_file, header=None, names=['Status'])
        nodes['Status'] = status['Status']
        
        # Calculate the distance and assign each hepatocyte to each zone using the 'q' parameter
        nodes['Distance'] = np.linalg.norm(nodes[['X', 'Y']].values - centre.values, axis=1)
        nodes['Zone'] = pd.qcut(nodes['Distance'], q=3, labels=['Zone 3', 'Zone 2', 'Zone 1'])      
        
        # Aggregate data by zone and status
        summary = nodes.groupby(['Zone', 'Status']).size().unstack(fill_value=0)
        summary = summary.reindex(columns=[0, 1, 2, 3], fill_value=0)
        summary.columns = ['Healthy', 'Damaged', 'Necrotic', 'Regenerated']
        summary_transposed = summary.transpose()  # J'ai tout transposé car une ligne parasite apparaissait à chaque fois

        ws = wb.create_sheet(title=f"Iteration {i}")
        # Write the transposed summary data into the Excel sheet
        for row in dataframe_to_rows(summary_transposed, index=True, header=True):
            ws.append(row)
    wb.save(output_excel_file_path)
    print ('Excel file generated')
    return

def generate_concentration_heatmap(substance, column_index, concentration_data_file, output_file, show_min_max, global_min, global_max):
   
    # Path to the file containing spatial coordinates of hepatocytes
    path_to_coordinates = 'Field_Nodes.txt'

    # Read the concentration data file, selecting the specified column for concentration values
    df_concentrations = pd.read_csv(concentration_data_file, sep=',', header=None, usecols=[0, column_index])
    df_concentrations.columns = ['HepatocyteID', f'{substance}_Concentration']
    df_concentrations['HepatocyteID'] = df_concentrations['HepatocyteID'].astype(int)

    # Read the coordinates data file
    df_coordinates = pd.read_csv(path_to_coordinates, sep=',', header=None, names=['HepatocyteID', 'X', 'Y'])
    df_coordinates['HepatocyteID'] = df_coordinates['HepatocyteID'].astype(int)

    # Merge the concentration data with the coordinates data
    df_merged = pd.merge(df_concentrations, df_coordinates, on='HepatocyteID')
    
    # Prepare data for interpolation
    resolution = 1000
    x, y, z = df_merged['X'].values, df_merged['Y'].values, df_merged[f'{substance}_Concentration'].values
    xi, yi = np.linspace(x.min(), x.max(), resolution), np.linspace(y.min(), y.max(), resolution)
    xi, yi = np.meshgrid(xi, yi)
    zi = griddata((x, y), z, (xi, yi), method='linear')

    # Generate the heatmap
    fig, ax = plt.subplots(figsize=(12, 8))
    scale = np.linspace(global_min, global_max, 1000, endpoint=True)
    if show_min_max:
        contour = ax.contourf(xi, yi, zi, scale, vmin=global_min, vmax=global_max, extend='both')
        plt.colorbar(contour, label=f'{substance} Concentration')
    else:
        contour = ax.contourf(xi, yi, zi, levels=100, cmap='viridis')
        plt.colorbar(contour, label=f'{substance} Concentration', ax=ax)

    # Add a red circle to represent the central vein if required
    center_x, center_y = np.mean(x), np.mean(y)
    circle_radius = 0.00000005 * resolution  
    small_circle_radius = 0.00000002 * resolution
    
    circle = plt.Circle((center_x, center_y), circle_radius, facecolor='blue', edgecolor='blue', fill=True, linewidth=2)
    ax.add_patch(circle)
    portal_vein1 = plt.Circle((np.mean(x), np.max(y)), small_circle_radius, facecolor='red', edgecolor='red', fill=True)
    ax.add_patch(portal_vein1)
    portal_vein2 = plt.Circle((np.mean(x), -np.max(y)), small_circle_radius, facecolor='red', edgecolor='red', fill=True)
    ax.add_patch(portal_vein2)
    portal_vein3 = plt.Circle((-np.max(x), np.max(y)/2), small_circle_radius, facecolor='red', edgecolor='red', fill=True)
    ax.add_patch(portal_vein3)
    portal_vein4 = plt.Circle((-np.max(x), -np.max(y)/2), small_circle_radius, facecolor='red', edgecolor='red', fill=True)
    ax.add_patch(portal_vein4)
    portal_vein5 = plt.Circle((np.max(x), -np.max(y)/2), small_circle_radius, facecolor='red', edgecolor='red', fill=True)
    ax.add_patch(portal_vein5)
    portal_vein6 = plt.Circle((np.max(x), np.max(y)/2), small_circle_radius, facecolor='red', edgecolor='red', fill=True)
    ax.add_patch(portal_vein6)


    # Optionally display global min and max concentration values
    if show_min_max:
        plt.text(0.95, 0.01, f'Global Min: {global_min}\nGlobal Max: {global_max}',
                 verticalalignment='bottom', horizontalalignment='right',
                 transform=ax.transAxes,
                 color='black', fontsize=10)

    plt.axis('off')
    # Save the generated heatmap image
    plt.savefig(output_file, dpi=300)
    plt.close()
    return

def generate_heatmaps_for_multiple_files(substance, column_index, concentration_files_base_path, output_files_base_path, num_files, show_min_max=True):
    last_iteration = num_files-1
    global_min, global_max = float('inf'), float('-inf')

    # Determine global min and max concentrations across all files
    for i in range(0, last_iteration + 1, 6):
        concentration_data_path = f"{concentration_files_base_path}-{i}.txt"
        df_concentrations = pd.read_csv(concentration_data_path, sep=',', header=None, usecols=[column_index])
        current_min = df_concentrations.iloc[column_index].min()
        current_max = df_concentrations.iloc[column_index].max()
        global_min, global_max = min(global_min, current_min), max(global_max, current_max)
    
    # Generate and save a heatmap for each file
    for i in range(0, last_iteration + 1, 6):
        concentration_data_path = f"{concentration_files_base_path}-{i}.txt"
        output_file = f"{output_files_base_path}_heatmap_{substance}_{i}.png"
        generate_concentration_heatmap(substance, column_index, concentration_data_path, output_file, show_min_max, global_min, global_max)
        print(f"Generated image for {output_file} / {last_iteration}")
    concentration_data_path = f"{concentration_files_base_path}-{last_iteration}.txt"
    output_file = f"{output_files_base_path}_heatmap_{substance}_{last_iteration}.png"
    generate_concentration_heatmap(substance, column_index, concentration_data_path, output_file, show_min_max, global_min, global_max)
    print(f"Generated image for {output_file} / {last_iteration}")
    return

def plot_combined_layer_with_legend_and_save(nodes_file, status_file, output_file):
    # Load the data
    nodes = pd.read_csv(nodes_file, sep=',', header=None, names=['Hepatocyte', 'X', 'Y'])
    status = pd.read_csv(status_file, header=None, names=['Status'])

    # Merge the datasets
    data = pd.concat([nodes, status], axis=1)

    # Increase resolution for better quality
    resolution = 1000  # Higher resolution for better image quality

    # Define a high-resolution grid to interpolate
    grid_x, grid_y = np.mgrid[min(data['X']):max(data['X']):resolution*1j, min(data['Y']):max(data['Y']):resolution*1j]

    # Define custom color maps for each status
    status_colors = ['green', 'orange', 'black', 'gray']
    status_maps = [LinearSegmentedColormap.from_list(f'status_{i}_cmap', ['white', color], N=2) for i, color in enumerate(status_colors)]

    # Initialize an array for each separate image and for the combined image
    separate_images = []
    combined_image = np.full((resolution, resolution, 3), 255)  # Start with a white image

    # Iterate over each status and create an individual layer
    for status_value, cmap in zip(range(4), status_maps):
        # Interpolate the data for the current status
        grid_z = griddata(data[['X', 'Y']], data['Status'] == status_value, (grid_x, grid_y), method='cubic', fill_value=0)

        # Plot the current status layer
        im = plt.imshow(grid_z.T, extent=(min(data['X']), max(data['X']), min(data['Y']), max(data['Y'])),
                        origin='lower', cmap=cmap, alpha=1.0)
        plt.clf()  # Clear the plot to avoid showing it

        # Extract the color image
        im.set_clim(0, 1)
        color_image = im.to_rgba(grid_z.T, bytes=True, norm=True)[:, :, :3]  # Remove the alpha channel

        # Add the current layer to the separate images list
        separate_images.append(color_image)

    # Determine the dominant color for each pixel in the combined image
    for i in range(resolution):
        for j in range(resolution):
            for img in separate_images:
                if not np.all(img[i, j] == 255):  # If the pixel is not white
                    combined_image[i, j] = img[i, j]
                    break

    # Plot the combined image with legend and circle
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.imshow(combined_image)
    ax.set_title('Hepatocyte Status Visualization')

   
    circle_center = (resolution/2, resolution/2)  
    circle_radius = 0.05 * resolution  
    small_circle_radius = 0.02 * resolution  

    central_vein = Circle(circle_center, circle_radius, facecolor='blue', edgecolor='blue', fill=True)
    ax.add_patch(central_vein)
    
    portal_vein1 = Circle((resolution/2, resolution), small_circle_radius, facecolor='red', edgecolor='red', fill=True)
    ax.add_patch(portal_vein1)
    portal_vein2 = Circle((resolution/2, 0), small_circle_radius, facecolor='red', edgecolor='red', fill=True)
    ax.add_patch(portal_vein2)
    portal_vein3 = Circle((0, resolution/4), small_circle_radius, facecolor='red', edgecolor='red', fill=True)
    ax.add_patch(portal_vein3)
    portal_vein4 = Circle((0, 3*resolution/4), small_circle_radius, facecolor='red', edgecolor='red', fill=True)
    ax.add_patch(portal_vein4)
    portal_vein5 = Circle((resolution, resolution/4), small_circle_radius, facecolor='red', edgecolor='red', fill=True)
    ax.add_patch(portal_vein5)
    portal_vein6 = Circle((resolution, 3*resolution/4), small_circle_radius, facecolor='red', edgecolor='red', fill=True)
    ax.add_patch(portal_vein6)


    # Create a custom legend
    custom_labels = ['Healthy (Green)', 'Damaged (Orange)', 'Necrotic (Black)', 'Regenerated (Gray)']
    custom_colors = [status_maps[i](1) for i in range(4)]
    patches = [plt.plot([], [], marker="o", ms=10, ls="", mec=None, color=color, 
                        label="{:s}".format(label))[0] for color, label in zip(custom_colors, custom_labels)]
    plt.axis('off')
    ax.legend(handles=patches, loc='upper right')
    
    # ax.set_xlabel("X Coordinate")
    # ax.set_ylabel("Y Coordinate")

    plt.savefig(output_file, dpi=300, bbox_inches='tight')  # Save with high DPI for
    plt.close()
    return
    
def generate_all_status_images(nodes_file, status_files_base_path, output_files_base_path, num_files):
    last_iteration = num_files-1
    for i in range(0, num_files, 6):
        # Construct file paths for status files and output images
        status_file = f"{status_files_base_path}-{i}.txt"
        output_file = f"{output_files_base_path}_{i}.png"

        # Generate and save the imag
        plot_combined_layer_with_legend_and_save(nodes_file, status_file, output_file)
        print(f"Generated image for {output_file} / {last_iteration}")
    status_file = f"{status_files_base_path}-{last_iteration}.txt"
    output_file = f"{output_files_base_path}_{last_iteration}.png"
    # Generate and save the imag
    plot_combined_layer_with_legend_and_save(nodes_file, status_file, output_file)
    print(f"Generated image for {output_file} / {last_iteration}")
    return



def OutputValuesStatus (iterations, fileNamePrefixe):
    '''
    Create plot of the % of each type of cell and saves it in a file

    Parameters
    ----------
    iterations : int
        Number of itteration.
    fileNamePrefixe : str

    Returns
    -------
    None.

    '''
    
    hepatocytesStatus_array  = np.zeros((iterations, 4))
    
    for i in range (iterations):
        
        hepatocytesStatusFile = open(fileNamePrefixe +'/_hepatocytes_status_Hepatocytes-%d.txt' %(i) , 'r+')
        hepatocyteStatus = hepatocytesStatusFile.readlines()
        
        numTotal = len(hepatocyteStatus)
        numSain = 0
        numDamaged = 0
        numDead = 0
        numRegen = 0
        
        for status in hepatocyteStatus:
            if status == '0.0\n':
                numSain = numSain + 1
            elif status =='1.0\n':
                numDamaged = numDamaged + 1
            elif status =='2.0\n':
                numDead = numDead + 1
            else:
                numRegen = numRegen + 1
        
        percentSain = 100 * numSain / numTotal
        percentDamaged = 100 * numDamaged / numTotal
        percentDead = 100 * numDead / numTotal
        percentRegen = 100 * numRegen / numTotal
        
        hepatocytesStatus_array[i][0] = percentSain
        hepatocytesStatus_array[i][1] = percentDamaged
        hepatocytesStatus_array[i][2] = percentDead
        hepatocytesStatus_array[i][3] = percentRegen
        
        # print('for iteration number ' + str(i))
        # print('% of sain = ' + str(percentSain))
        # print('% of damaged = ' + str(percentDamaged))
        # print('% of dead = ' + str(percentDead))
        
    hepatocytesSain  = np.zeros(iterations)
    hepatocytesDam  = np.zeros(iterations)
    hepatocytesDead  = np.zeros(iterations)
    hepatocytesRegen  = np.zeros(iterations)
    for i in range (len(hepatocytesStatus_array)):
        hepatocytesSain[i] = hepatocytesStatus_array[i][0]
        hepatocytesDam[i] = hepatocytesStatus_array[i][1]
        hepatocytesDead[i] = hepatocytesStatus_array[i][2]
        hepatocytesRegen[i] = hepatocytesStatus_array[i][3]
    fig, ax = plt.subplots()
    ax.plot(hepatocytesSain, 'g', label='Healthy')
    ax.plot(hepatocytesDam, 'r', label='Damaged')
    ax.plot(hepatocytesDead, 'k', label='Necrosed')
    ax.plot(hepatocytesRegen, 'm', label='Regenerated')
    plt.xlabel("time (*20min)")
    plt.ylabel("% of cells")
    plt.suptitle("% of cell at each status in the lobule")
    plt.legend()
    # plt.show()
    plt.savefig(fileNamePrefixe + "/Plot_status.png")
    plt.close()
    print(f"Generated status plot for {fileNamePrefixe}")
    return 


def OutputValuesCells (iterations, fileNamePrefixe, valueType, nodes):
    '''
    Create plot of the evolution of values on different nodes and saves it in a file

    Parameters
    ----------
    iterations : int
        Number of itteration.
    fileNamePrefixe : str
    
    valueType : int
        1 = APAP ; 2 = NAPQI ; 3 = GSH ; 4 = Healthy ; 5 = Damaged ; 6 = Necrosed ; 7 = Regen ; 8 = AST ; 9 = ALT.
    nodes : Array
        Array of the node number to look at.

    Returns
    -------
    hepatocytesValues_array : array
        Array of values with every itteration value for every nodes.

    '''
    # hepatocytesValues_array  = numpy.zeros((iterations, len(nodes)))
    fig, ax = plt.subplots()
    variable_Equivalent = ['quantity of APAP', 'quantity of NAPQI', 'quantity of GSH', 'chances of healthy', 'chances of damaged', 'chances of necrosed', 'chances of regeneration', 'AST concentration', 'ALT concentration']
    n = 0
    for node in nodes:
        hepatocytesValues_array  = np.zeros(iterations)
        for i in range (iterations):
            hepatocytesValuesFile = open(fileNamePrefixe +'/_Used-Field_Hepatocytes-%d.txt' %(i) , 'r+')
            hepatocyteCellValues = hepatocytesValuesFile.readlines()
            valuesNode = hepatocyteCellValues[node-1].split(',')
            # hepatocytesValues_array[i][n] = valuesNode[valueType]
            hepatocytesValues_array[i] = valuesNode[valueType]
            hepatocytesValuesFile.close()
        n=n+1
        ax.plot(hepatocytesValues_array, label='Node ' + str(node))
    
    plt.xlabel("time (*20min)")
    plt.ylabel(variable_Equivalent[valueType-1])
    plt.suptitle("Evolution of the " + variable_Equivalent[valueType-1])
    plt.legend()
    # plt.show()
    plt.savefig(fileNamePrefixe + '/Plot_' + variable_Equivalent[valueType-1] +".png")
    plt.close()
    print(f"Generated value {valueType} plot for {fileNamePrefixe}")
    return 

def OutputValuesASTALT (iterations, fileNamePrefixe):
    
    nodeField = open('Field_Nodes.txt' , 'r+')
    resultsFile = open('ResultsASTALTValues.txt', 'a+')
    nodes = nodeField.readlines()
    hepatocytesASTValues_array  = np.zeros(iterations)
    hepatocytesALTValues_array  = np.zeros(iterations)

    for i in range (iterations):
        hepatocytesValuesFile = open(fileNamePrefixe +'/_Used-Field_Hepatocytes-%d.txt' %(i) , 'r+')
        hepatocyteCellValues = hepatocytesValuesFile.readlines()
        for n in range(len(nodes)):
            valuesNode = hepatocyteCellValues[n-1].split(',')
            hepatocytesASTValues_array[i] = hepatocytesASTValues_array[i] + float(valuesNode[8])
            hepatocytesALTValues_array[i] = hepatocytesALTValues_array[i] + float(valuesNode[9])
        hepatocytesASTValues_array[i] = hepatocytesASTValues_array[i] / (len(nodes))
        hepatocytesALTValues_array[i] = hepatocytesALTValues_array[i] / (len(nodes))
        hepatocytesValuesFile.close()
        if i == iterations-1:
            resultsFile.write(str(fileNamePrefixe) + '\t' + str(hepatocytesASTValues_array[i]) + '\t' + str(hepatocytesALTValues_array[i]) + '\n')
        
    time = np.linspace(0, 20*iterations, iterations)
    
    fig, axs = plt.subplots(2, 1,figsize=(20, 10))
    fig.suptitle('Evolution of the concentration of AST and ALT in the lobule')
    axs[0].plot(time, hepatocytesASTValues_array)
    axs[0].set_title('Concentration of AST (IU/L)')
    axs[1].plot(time, hepatocytesALTValues_array)
    axs[1].set_title('Concentration of ALT (IU/L)')
    for ax in axs.flat:
        ax.set(xlabel='Time (s)', ylabel='')

    fig.tight_layout()
    fig.savefig(fileNamePrefixe + "/Plot_AST-ALT.png")
    plt.close()
    nodeField.close()
    resultsFile.close()
    print(f"Generated AST/ALT plot for {fileNamePrefixe}")
    return


## TO USE AS A STANDALONE AFTER A SIMULATION HAS BEEN RUN AND WANTING TO OUTPUT ALL RESULTS
def Print_All_Results(input_file, path_prefix):
    # Example usage with your file paths
    nodes_file_path = 'Field_Nodes.txt'
    
    csvfile = pd.read_excel(input_file).values
    # for c in range(2, 3):
    for c in range(2, len(csvfile[0])):
        patientNum = str(int(csvfile[0][c]))
        totalIterations = int(csvfile[1][c] * 72)
        
        folder_base_path = path_prefix + str(patientNum)
        status_files_base_path = folder_base_path + '/_hepatocytes_status_Hepatocytes'
        concentration_files_base_path = folder_base_path + '/_Used-Field_Hepatocytes'
        output_files_base_path = folder_base_path + '/output_image'
        output_excel_file_path = folder_base_path + '/hepatocyte_zone_stats.xlsx'
    
        generate_all_status_images(nodes_file_path, status_files_base_path, output_files_base_path, totalIterations )
        generate_heatmaps_for_multiple_files('APAP', 1, concentration_files_base_path, output_files_base_path, totalIterations, True)
        generate_heatmaps_for_multiple_files('NAPQI', 2, concentration_files_base_path, output_files_base_path, totalIterations, True)
        generate_heatmaps_for_multiple_files('GSH', 3, concentration_files_base_path, output_files_base_path, totalIterations, True)
        generate_heatmaps_for_multiple_files('AST', 8, concentration_files_base_path, output_files_base_path, totalIterations, True)
        generate_heatmaps_for_multiple_files('ALT', 9, concentration_files_base_path, output_files_base_path, totalIterations, True)
        generate_excel(nodes_file_path, status_files_base_path, output_excel_file_path, totalIterations)
        OutputValuesCells(totalIterations, folder_base_path, 1, [2494, 371, 1])
        OutputValuesCells(totalIterations, folder_base_path, 2, [2494, 371, 1])
        OutputValuesCells(totalIterations, folder_base_path, 3, [2494, 371, 1])
        OutputValuesCells(totalIterations, folder_base_path, 4, [2494, 371, 1])
        OutputValuesCells(totalIterations, folder_base_path, 5, [2494, 371, 1])
        OutputValuesCells(totalIterations, folder_base_path, 6, [2494, 371, 1])
        OutputValuesCells(totalIterations, folder_base_path, 7, [2494, 371, 1])
        
        OutputValuesStatus(totalIterations, folder_base_path)

        OutputValuesASTALT(totalIterations, folder_base_path)
    return


## TO USE WHEN RUNNING A SIMULATION AND WANTING TO OUTPUT ALL RESULTS
def Print_Results(folder_base_path, totalIterations):
    # Example usage with your file paths
    nodes_file_path = 'Field_Nodes.txt'
    
    status_files_base_path = folder_base_path + '/_hepatocytes_status_Hepatocytes'
    concentration_files_base_path = folder_base_path + '/_Used-Field_Hepatocytes'
    output_files_base_path = folder_base_path + '/output_image'
    output_excel_file_path = folder_base_path + '/hepatocyte_zone_stats.xlsx'

    generate_all_status_images(nodes_file_path, status_files_base_path, output_files_base_path, totalIterations )
    
    generate_heatmaps_for_multiple_files('APAP', 1, concentration_files_base_path, output_files_base_path, totalIterations, True)
    generate_heatmaps_for_multiple_files('NAPQI', 2, concentration_files_base_path, output_files_base_path, totalIterations, True)
    generate_heatmaps_for_multiple_files('GSH', 3, concentration_files_base_path, output_files_base_path, totalIterations, True)
    generate_heatmaps_for_multiple_files('AST', 8, concentration_files_base_path, output_files_base_path, totalIterations, True)
    generate_heatmaps_for_multiple_files('ALT', 9, concentration_files_base_path, output_files_base_path, totalIterations, True)
    
    generate_excel(nodes_file_path, status_files_base_path, output_excel_file_path, totalIterations)
    
    OutputValuesCells(totalIterations, folder_base_path, 1, [2494, 371, 1])
    OutputValuesCells(totalIterations, folder_base_path, 2, [2494, 371, 1])
    OutputValuesCells(totalIterations, folder_base_path, 3, [2494, 371, 1])
    OutputValuesCells(totalIterations, folder_base_path, 4, [2494, 371, 1])
    OutputValuesCells(totalIterations, folder_base_path, 5, [2494, 371, 1])
    OutputValuesCells(totalIterations, folder_base_path, 6, [2494, 371, 1])
    OutputValuesCells(totalIterations, folder_base_path, 7, [2494, 371, 1])
    
    OutputValuesStatus(totalIterations, folder_base_path)

    OutputValuesASTALT(totalIterations, folder_base_path)
    return

# Print_All_Results('InputValues.xlsx', 'Results-')