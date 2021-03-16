#   ---RAO NOZZLE DESIGNER---
#   METU ROCKET SOCIETY, 2021

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tkinter import *

root = Tk()
root.title("RAO NOZZLE DESIGNER")
root.iconbitmap('mrs.ico')

throat_radius = 0
exit_radius = 0
chamber_radius = 0
expansion_ratio = 0
delta_n = 0
theta_FC = 0

#Input fields, with default values
throat_radius_label = Label(root, text="Throat Radius (mm): ").grid(row=2, column=0)
throat_radius_input = Entry(root)
throat_radius_input.grid(row=2, column=1)
throat_radius_input.insert(0, "24.")

exit_radius_label = Label(root, text="Exit Radius (mm): ").grid(row=3, column=0)
exit_radius_input = Entry(root)
exit_radius_input.grid(row=3, column=1)
exit_radius_input.insert(0, "41.5")

chamber_radius_label = Label(root, text="Chamber Radius (mm): ").grid(row=4, column=0)
chamber_radius_input = Entry(root)
chamber_radius_input.grid(row=4, column=1)
chamber_radius_input.insert(0, "48.2")

expansion_ratio_label = Label(root, text="Expansion Ratio: ").grid(row=5, column=0)
expansion_ratio_input = Entry(root)
expansion_ratio_input.grid(row=5, column=1)
expansion_ratio_input.insert(0, "2.96")

delta_n_label = Label(root, text="Delta n: ").grid(row=6, column=0)
delta_n_input = Entry(root)
delta_n_input.grid(row=6, column=1)
delta_n_input.insert(0, "22.")

theta_FC_label = Label(root, text="Theta FC: ").grid(row=7, column=0)
theta_FC_input = Entry(root)
theta_FC_input.grid(row=7, column=1)
theta_FC_input.insert(0, "-90.")

first_curve_dict = {}
second_curve_dict = {}
parabolic_curve_dict = {}

#######################
def nozzle_length(throat_radius, expansion_ratio):
    Ln = 0.80*(math.sqrt(expansion_ratio)-1)*throat_radius/math.tan(math.radians(15))
    return Ln

def firstcurve_xfc (theta_FC, throat_radius):
    x_fc = float( math.cos(math.radians(theta_FC))*1.5*throat_radius)
    return x_fc

def firstcurve_yfc (theta_FC, throat_radius):
    y_fc = float( math.sin(math.radians(theta_FC))*1.5*throat_radius + (1.5*throat_radius+throat_radius))
    return y_fc

def secondcurve_xsc (theta_SC, throat_radius):
    x_sc = float( math.cos(math.radians(theta_SC))*0.382*throat_radius)
    return x_sc

def secondcurve_ysc (theta_SC, throat_radius):
    y_sc = float(math.sin(math.radians(theta_SC))*0.382*throat_radius + (0.382*throat_radius+throat_radius))
    return y_sc

def paraboliccurvematrix (Y_exit,X_exit, delta_n, X_SC, Y_SC):
    third_invariant = 1/(math.tan(math.radians(delta_n)))
    arr = np.array([[Y_SC**2, Y_SC, 1.0], [Y_exit**2, Y_exit, 1.0], [Y_SC*2, 1.0, 0.0]])
    para = np.linalg.inv(arr)
    arr2 = np.array([X_SC, X_exit, third_invariant])
    i_matrix = para.dot(arr2)
    return i_matrix

def exit_angle(Ln, throat_radius,exit_radius):
    dy = exit_radius-throat_radius
    delta_e = float (math.degrees(math.atan(dy/Ln)))
    return delta_e
#######################

def calculateGeometry(excelFilename): #INITIALIZE NOZZLE GEOMETRY CALCULATION
    print("Calculating...")
    try:
        throat_radius = float(throat_radius_input.get())
        exit_radius = float(exit_radius_input.get())
        chamber_radius = float(chamber_radius_input.get())
        expansion_ratio = float(expansion_ratio_input.get())
        delta_n = float(delta_n_input.get())
        theta_FC = float(theta_FC_input.get())
    except:
        print("Please enter floating point values in all inputs.")

    # Make sure we clean up dicts before running the next calculation
    first_curve_dict = {}
    second_curve_dict = {}
    parabolic_curve_dict = {}

    # First curve calculations

    theta_initial = math.degrees(math.asin((chamber_radius - (1.5 * throat_radius + throat_radius))/(1.5*throat_radius)))
    fc_steps = (-90-theta_initial)/20

    while True:
        x_fc = firstcurve_xfc(theta_FC = theta_FC, throat_radius = throat_radius)
        y_fc = firstcurve_yfc(theta_FC = theta_FC, throat_radius = throat_radius)
        first_curve_dict.update({ x_fc : y_fc })
        theta_FC = theta_FC + fc_steps

        if y_fc >= chamber_radius :
            theta_FC = theta_FC - 2*fc_steps
            break

    # Second curve calculations

    theta_SC = math.degrees(math.asin((throat_radius - (0.382 * throat_radius + throat_radius))/(0.382*throat_radius)))
    Sc_steps = (delta_n)/22
    while True:

        x_sc = secondcurve_xsc(theta_SC = theta_SC, throat_radius = throat_radius)
        y_sc = secondcurve_ysc(theta_SC = theta_SC, throat_radius = throat_radius)
        second_curve_dict.update({ x_sc : y_sc })
        theta_SC = theta_SC + Sc_steps

        if theta_SC >= (delta_n - 90) :
            theta_SC = theta_SC - 2*Sc_steps

            break

    # Exit angle

    Ln = nozzle_length(throat_radius, expansion_ratio)
    delta_e = exit_angle(Ln, throat_radius, exit_radius)

    # Parabolic curve

    second_curve_x = max(second_curve_dict.keys())
    second_curve_y = max(second_curve_dict.values())


    i_matrix = paraboliccurvematrix(exit_radius,Ln,delta_n,second_curve_x, second_curve_y )

    a = i_matrix[0]
    b = i_matrix[1]
    c = i_matrix[2]
    parabolic_x = list()
    parabolic_y = list()
    parabolic_step=(exit_radius-second_curve_y)/21
    while True:
        px = a*second_curve_y**2+b*second_curve_y+c
        parabolic_x.append(px)
        parabolic_y.append(second_curve_y)
        second_curve_y = second_curve_y + parabolic_step
        if second_curve_y > exit_radius:
            break

    if not excelFilename == "" or excelFilename == None:
        if len(excelFilename) > 5 and excelFilename[-5] == ".":
            exportFile = excelFilename
        else:
            exportFile = excelFilename + ".xlsx"
            
        rao_nozzle = {'Xfc': list(first_curve_dict.keys()),'Yfc': list(first_curve_dict.values()), 'Xsc': list(second_curve_dict.keys()), 'Ysc': list(second_curve_dict.values()), 'Pxc': parabolic_x, 'Pyc': parabolic_y}
        rao_nozzle_df = pd.DataFrame(rao_nozzle)

        #rao_nozzle_df.to_excel("rao_nozzle_geometry.xlsx", sheet_name= "coordinates",)
        rao_nozzle_df.to_excel(exportFile, sheet_name= "coordinates",)

        inputSaveFile = exportFile[0:-5] + ".txt"
        result_file = open(inputSaveFile, "w")
        result_file.write("INPUTS")
        result_file.write(" throat radius: ")
        result_file.write( str(throat_radius))
        result_file.write(" chamber radius: ")
        result_file.write(str(chamber_radius))
        result_file.write(" exit radiues: ")
        result_file.write(str(exit_radius))
        result_file.write(" throat angle: ")
        result_file.write(str(delta_n))
        result_file.write(" exit angle: ")
        result_file.write(str(delta_e))
        result_file.close()

        print("Successfully saved geometry to " + exportFile)
        print("Inputs saved in " + inputSaveFile)

    else:
        print("Skipping export...")

    exit_angle_label = Label(root, text="Exit angle: " + str(delta_e))
    Label(root, text="").grid(row=13, column=0) #padding
    exit_angle_label.grid(row=14, column=0)

    print("Done.")

    #plt.clf()
    plt.figure()
    plt.title("Nozzle Geometry")
    plt.xlabel('X Axis (mm)')
    plt.xlim(-80,80)
    plt.ylabel('Y Axis (mm)')
    plt.ylim(0,60)

    plt.plot(list(first_curve_dict.keys()),list(first_curve_dict.values()), "b")
    plt.plot(list(second_curve_dict.keys()),list(second_curve_dict.values()),"r")
    plt.plot(parabolic_x,parabolic_y,"g")
    plt.show()
    

#Click button to calculate nozzle (and optionally export to Excel file)

Label(root, text="").grid(row=8, column=0) #padding, just for aesthetic purposes

excel_filename_label = Label(root, text="Geometry Export Filename (leave blank to skip export): ").grid(row=9, column=0, columnspan=2)
excel_filename = Entry(root, width=50)
excel_filename.grid(row=10, column=0, columnspan=2)
excel_filename.insert(0, "rao_nozzle_geometry.xlsx")

Label(root, text="").grid(row=11, column=0) #padding

runButton = Button(text="Calculate Nozzle Geometry", command=lambda:calculateGeometry(str(excel_filename.get())))
runButton.grid(row=12, column=0)

root.mainloop()
