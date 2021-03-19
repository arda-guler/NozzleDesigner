#   ---RAO NOZZLE DESIGNER---
#   METU ROCKET SOCIETY, 2021

from dearpygui.core import *
from dearpygui.simple import *
import math
import numpy as np
import pandas as pd

#set initial window configuration (purely cosmetic)
set_main_window_size(1000, 650)
set_main_window_title("Rao Nozzle Designer | MRS")
set_theme("Dark")

def importFile():
    try:
        import_filepath = get_value("load_path_field")
        log_info("Importing inputs from " + import_filepath, logger="Logs")
        import_file = open(import_filepath, "r")
    except:
        log_error("Import failed. Check filepath.", logger="Logs")
        return

    try:
        import_lines = import_file.readlines()
        set_value(name="throat_radius_field", value=import_lines[2][15:-1])
        set_value(name="exit_radius_field", value=import_lines[3][13:-1])
        set_value(name="chamber_radius_field", value=import_lines[4][16:-1])
        set_value(name="expansion_ratio_field", value=import_lines[5][17:-1])
        set_value(name="delta_n_field", value=import_lines[6][24:-1])
        set_value(name="theta_FC_field", value=import_lines[7][10:-1])
    except:
        log_error("Import failed. Check file formatting.", logger="Logs")
        return
    
    log_info("Import successful.", logger="Logs")
  

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
    log_info(message = "Computing nozzle geometry...", logger = "Logs")
    try:
        throat_radius = float(get_value("throat_radius_field"))
        exit_radius = float(get_value("exit_radius_field"))
        chamber_radius = float(get_value("chamber_radius_field"))
        expansion_ratio = float(get_value("expansion_ratio_field"))
        delta_n = float(get_value("delta_n_field"))
        theta_FC = float(get_value("theta_FC_field"))

        excelFilename = get_value("filename_field")
    except:
        log_error("Input error. Make sure all design parameters are float values and export filename doesn't include illegal characters.", logger = "Logs")

    log_info("Inputs:\n" +
             "Throat Radius: " + str(throat_radius) + " mm\n"
             "Exit Radius: " + str(exit_radius) + " mm\n"
             "Chamber Radius: " + str(chamber_radius) + " mm\n"
             "Expansion Ratio: " + str(expansion_ratio) + "\n"
             "Delta n: " + str(delta_n) + "\n"
             "Theta FC: " + str(theta_FC), logger = "Logs")

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


    i_matrix = paraboliccurvematrix(exit_radius, Ln,delta_n, second_curve_x, second_curve_y)

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

    # File output
    if not excelFilename == "" or excelFilename == None:
        log_info("Attempting export...", logger = "Logs")
        if len(excelFilename) > 5 and excelFilename[-5] == ".":
            exportFile = excelFilename
        else:
            exportFile = excelFilename + ".xlsx"

        try:           
            rao_nozzle_fc = {'X_First_Curve': list(first_curve_dict.keys())[::-1],'Y_First_Curve': list(first_curve_dict.values())[::-1]}
            # (first curve is inverted for excel export so that all x values are in increasing order)
            rao_nozzle_sc = {'X_Second_Curve': list(second_curve_dict.keys()), 'Y_Second_Curve': list(second_curve_dict.values())}
            rao_nozzle_pc = {'X_Parabolic_Curve': parabolic_x, 'Y_Parabolic_Curve': parabolic_y}

            rao_nozzle_fc_df = pd.DataFrame(rao_nozzle_fc)
            rao_nozzle_sc_df = pd.DataFrame(rao_nozzle_sc)
            rao_nozzle_pc_df = pd.DataFrame(rao_nozzle_pc)

            with pd.ExcelWriter(exportFile) as writer:
                rao_nozzle_fc_df.to_excel(writer, sheet_name = 'coordinates', startcol=0)
                rao_nozzle_sc_df.to_excel(writer, sheet_name = 'coordinates', startcol=3)
                rao_nozzle_pc_df.to_excel(writer, sheet_name = 'coordinates', startcol=6)
  
            log_info("Successfully saved geometry to " + exportFile, logger = "Logs")
            
        except:
            log_error("Excel export failed.", logger = "Logs")

        try:
            inputSaveFile = exportFile[0:-5] + ".txt"
            result_file = open(inputSaveFile, "w")
            result_file.write("INPUTS\n\n")
            result_file.write("Throat radius: ")
            result_file.write(str(throat_radius)+"\n")
            result_file.write("Exit radius: ")
            result_file.write(str(exit_radius)+"\n")
            result_file.write("Chamber radius: ")
            result_file.write(str(chamber_radius)+"\n")
            result_file.write("Expansion ratio: ")
            result_file.write(str(expansion_ratio)+"\n")
            result_file.write("Throat angle (delta n): ")
            result_file.write(str(delta_n)+"\n")
            result_file.write("Theta FC: ")
            result_file.write(str(get_value("theta_FC_field"))+"\n")
            result_file.write("\nOUTPUTS\n\n")
            result_file.write("Exit angle: ")
            result_file.write(str(delta_e)+"\n")
            result_file.close()
            log_info("Inputs saved in " + inputSaveFile, logger = "Logs")
        except:
            log_error("TXT export failed.", logger = "Logs")  

    else:
        log_info("Skipping export...", logger = "Logs")

    log_info("Done.", logger = "Logs")

    #add data to plot
    add_line_series(name="Upstream Curve" , plot="nozzle_geometry",x=list(first_curve_dict.keys()), y=list(first_curve_dict.values()))
    add_line_series(name="Init. Exp. Curve" , plot="nozzle_geometry",x=list(second_curve_dict.keys()), y=list(second_curve_dict.values()))
    add_line_series(name="Parabolic Curve" , plot="nozzle_geometry",x=parabolic_x, y=parabolic_y)

#INPUTS WINDOW
with window("Input", width=450, height=350):   
    set_window_pos("Input", 10, 10)
    add_text("Enter design parameters in float values.")
    add_spacing(count=6)
    add_input_text(name = "throat_radius_field", label = "Throat Radius (mm)")
    add_input_text(name = "exit_radius_field", label = "Exit Radius (mm)")
    add_input_text(name = "chamber_radius_field", label = "Chamber Radius (mm)")
    add_input_text(name = "expansion_ratio_field", label = "Expansion Ratio")
    add_input_text(name = "delta_n_field", label = "Delta n (Throat Angle, deg)")
    add_input_text(name = "theta_FC_field", label = "Theta FC")
    add_spacing(count=6)
    add_input_text(name = "filename_field", label = "Export Filename", tip = "Leave blank to skip export. File extension is automatic.")
    add_spacing(count=6)
    add_button("Compute Nozzle Geometry", callback = computeNozzle)
    add_spacing(count=6)
    add_text("Alternatively, you can import a previously saved .txt file.")
    add_spacing(count=6)
    add_input_text(name = "load_path_field", label = "Import Filepath", tip = "If the file is in the same directory with the script, you don't need\n to write the full path.")
    add_spacing(count=6)
    add_button("Import File", callback = importFile)

#OUTPUTS WINDOW
with window("Output", width=500, height=400):
    set_window_pos("Output", 470, 10)
    add_input_text(name="exit_angle_output", label="Exit Angle (deg)", source="exit_angle", readonly=True, enabled=False)
    add_plot(name="nozzle_geometry", label="Nozzle Geometry",
             x_axis_name="X Axis (mm)", y_axis_name = "Y Axis (mm)", equal_aspects = True)

#LOG WINDOW
with window("Log", width=450, height=200):
    set_window_pos("Log", 10, 370)
    add_logger("Logs", log_level=0, autosize_x = True, autosize_y = True)

start_dearpygui()
