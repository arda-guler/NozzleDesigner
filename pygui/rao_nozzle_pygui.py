#   ---RAO NOZZLE DESIGNER---
#   METU ROCKET SOCIETY, 2021

from dearpygui.core import *
from dearpygui.simple import *
import math
import numpy as np
import pandas as pd

#set initial window configuration (purely cosmetic)
set_main_window_size(1000, 600)
set_main_window_title("MRS - RAO NOZZLE DESIGNER")
set_theme("Dark")

#init values at start for safety measures
throat_radius = 0
exit_radius = 0
chamber_radius = 0
expansion_ratio = 0
delta_n = 0
theta_FC = 0

first_curve_dict = {}
second_curve_dict = {}
parabolic_curve_dict = {}

# Calculation sub-functions

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

################################
#   THE ACTUAL CALCULATIONS
################################

def computeNozzle():
    print("Computing nozzle geometry...")
    try:
        throat_radius = float(get_value("throat_radius_field"))
        exit_radius = float(get_value("exit_radius_field"))
        chamber_radius = float(get_value("chamber_radius_field"))
        expansion_ratio = float(get_value("expansion_ratio_field"))
        delta_n = float(get_value("delta_n_field"))
        theta_FC = float(get_value("theta_FC_field"))

        excelFilename = get_value("filename_field")
    except:
        print("Input error.")

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

    # Update exit angle to display
    set_value(name = "exit_angle", value = delta_e)

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

    print("Done.")

    #add data to plot
    add_line_series(name="Upstream Curve" , plot="nozzle_geometry",x=list(first_curve_dict.keys()), y=list(first_curve_dict.values()))
    add_line_series(name="Init. Exp. Curve" , plot="nozzle_geometry",x=list(second_curve_dict.keys()), y=list(second_curve_dict.values()))
    add_line_series(name="Parabolic Curve" , plot="nozzle_geometry",x=parabolic_x, y=parabolic_y)

#INPUTS WINDOW
with window("Inputs", width=450, height=350):
    
    #print("Program running.") #ultimate debug strategy :)
    
    set_window_pos("Inputs", 10, 10)
    add_text("Enter floating point values in all fields.")
    add_spacing(count=6)
    add_input_text(name = "throat_radius_field", label = "Throat Radius (mm)")
    add_input_text(name = "exit_radius_field", label = "Exit Radius (mm)")
    add_input_text(name = "chamber_radius_field", label = "Chamber Radius (mm)")
    add_input_text(name = "expansion_ratio_field", label = "Expansion Ratio")
    add_input_text(name = "delta_n_field", label = "Delta n")
    add_input_text(name = "theta_FC_field", label = "Theta FC")
    add_spacing(count=6)
    add_input_text(name = "filename_field", label = "Export Filename")
    add_spacing(count=6)
    add_button("Compute Nozzle Geometry", callback = computeNozzle)

#OUTPUTS WINDOW
with window("Output", width=500, height=400):
    set_window_pos("Output", 470, 10)
    add_input_text(name="exit_angle_output", label="Exit Angle (deg)", source="exit_angle", readonly=True, enabled=False)
    add_plot(name="nozzle_geometry", label="Nozzle Geometry",
             x_axis_name="X Axis (mm)", y_axis_name = "Y Axis (mm)")

start_dearpygui()
