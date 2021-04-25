#   ---RAO NOZZLE DESIGNER---
#   METU ROCKET SOCIETY, 2021
#
#   Version 1.1.1

from dearpygui.core import *
from dearpygui.simple import *
import math
import numpy as np
import pandas as pd

#set initial window configuration (purely cosmetic)
set_main_window_size(1000, 700)
set_main_window_title("Rao Nozzle Designer | MRS")
set_theme("Dark")

calc_run_number = 0

#variables to save values of last run
#saving in another variable in case user makes changes to the input fields before clicking Export
last_throat_radius = None
last_exit_radius = None
last_chamber_radius = None
last_delta_n = None
last_theta_FC = None

last_first_curve = None
last_second_curve = None
last_parabolic_curve = None

def importFile():
    try:
        import_filepath = get_value("filepath_field")
        
        if not import_filepath[-4:] == ".txt":
            if import_filepath[-5:] == ".xlsx":
                log_warning("Exported .xlsx files don't contain input info. Trying " + import_filepath[:-5] + ".txt instead...", logger="Logs")
                import_filepath = import_filepath[:-5] + ".txt"
            else:
                import_filepath = import_filepath + ".txt"
            
        log_info("Importing inputs from " + import_filepath, logger="Logs")
        import_file = open(import_filepath, "r")
    except:
        log_error("Import failed. Check filepath.", logger="Logs")
        return

    try:
        import_lines = import_file.readlines()
        if not import_lines[0][18:-1] == "1.1.1":
            log_warning("Save file version does not match software version. Import might fail.", logger="Logs")
        
        set_value(name="throat_radius_field", value=import_lines[4][15:-1])
        set_value(name="exit_radius_field", value=import_lines[5][13:-1])
        set_value(name="chamber_radius_field", value=import_lines[6][16:-1])
        set_value(name="delta_n_field", value=import_lines[7][24:-1])
        set_value(name="theta_FC_field", value=import_lines[8][10:-1])
    except:
        log_error("Import failed. Check file formatting.", logger="Logs")
        return

    log_info("Import successful.", logger="Logs")
    computeNozzle()

def exportFile():
    if not calc_run_number > 0:
        log_error("Cannot export. Run the calculations first.", logger="Logs")
        return
    
    excelFilename = get_value("filepath_field")
    
    if not excelFilename == "" or excelFilename == None:
        log_info("Attempting export...", logger = "Logs")
        if len(excelFilename) > 5 and excelFilename[-5:] == ".xlsx":
            exportFile = excelFilename
        elif len(excelFilename) > 4 and excelFilename[-4:] == ".txt":
            exportFile = excelFilename[:-4] + ".xlsx"
        else:
            exportFile = excelFilename + ".xlsx"

        try:
            rao_nozzle_fc = {'X_First_Curve': last_first_curve[0][::-1],'Y_First_Curve': last_first_curve[1][::-1]}
            # (first curve is inverted for excel export so that all x values are in increasing order)
            rao_nozzle_sc = {'X_Second_Curve': last_second_curve[0], 'Y_Second_Curve': last_second_curve[1]}
            rao_nozzle_pc = {'X_Parabolic_Curve': last_parabolic_curve[0], 'Y_Parabolic_Curve': last_parabolic_curve[1]}

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
            result_file.write("Save file version 1.1.1\n\n")
            result_file.write("INPUTS\n\n")
            result_file.write("Throat radius: ")
            result_file.write(str(last_throat_radius)+"\n")
            result_file.write("Exit radius: ")
            result_file.write(str(last_exit_radius)+"\n")
            result_file.write("Chamber radius: ")
            result_file.write(str(last_chamber_radius)+"\n")
            result_file.write("Throat angle (delta n): ")
            result_file.write(str(last_delta_n)+"\n")
            result_file.write("Theta FC: ")
            result_file.write(str(last_theta_FC)+"\n")
            result_file.write("\nOUTPUTS\n\n")
            result_file.write("Exit angle: ")
            result_file.write(str(get_value("exit_angle"))+"\n")
            result_file.write("Nozzle expansion ratio: ")
            result_file.write(str(get_value("nozzle_expansion_ratio"))+"\n")
            result_file.write("Chamber contraction ratio: ")
            result_file.write(str(get_value("chamber_contraction_ratio"))+"\n")
            result_file.write("Geometry export file: " + exportFile + "\n")
            result_file.close()
            log_info("Inputs saved in " + inputSaveFile, logger = "Logs")
        except:
            log_error("TXT export failed.", logger = "Logs")  

    else:
        log_info("Skipping export...", logger = "Logs")

    log_info("Done.", logger = "Logs")

def computeNozzle():
    global calc_run_number
    calc_run_number += 1
    log_info(message = "Run [" + str(calc_run_number) + "]: Computing nozzle geometry...", logger = "Logs")
    
    try:
        throat_radius = float(get_value("throat_radius_field"))
        exit_radius = float(get_value("exit_radius_field"))
        chamber_radius = float(get_value("chamber_radius_field"))
        delta_n = float(get_value("delta_n_field"))
        theta_FC = float(get_value("theta_FC_field"))
    except:
        log_error("Input error. Make sure all design parameters are float values.", logger = "Logs")

    # Right away we can calculate expansion and contraction ratios
    expansion_ratio = (exit_radius**2)/(throat_radius**2)
    set_value(name="nozzle_expansion_ratio", value=expansion_ratio)
    last_expansion_ratio = expansion_ratio

    chamber_contraction_ratio = (chamber_radius**2)/(throat_radius**2)
    set_value(name="chamber_contraction_ratio", value=chamber_contraction_ratio)
    last_contraction_ratio = chamber_contraction_ratio

    # save these values in global scope, in case we want to export
    global last_throat_radius, last_exit_radius, last_chamber_radius, last_delta_n, last_theta_FC
    last_throat_radius = throat_radius
    last_exit_radius = exit_radius
    last_chamber_radius = chamber_radius
    last_delta_n = delta_n
    last_theta_FC = theta_FC

    log_info("Inputs:\n" +
             "Throat Radius: " + str(throat_radius) + " mm\n"
             "Exit Radius: " + str(exit_radius) + " mm\n"
             "Chamber Radius: " + str(chamber_radius) + " mm\n"
             "Delta n: " + str(delta_n) + "\n"
             "Theta FC: " + str(theta_FC), logger = "Logs")

    # Make sure we clean up dicts before running the next calculation
    
    first_curve_dict = {}
    second_curve_dict = {}
    parabolic_curve_dict = {}

    ################################
    #   THE ACTUAL CALCULATIONS
    ################################

    # Calculation sub-functions

    def clamp(num, min_value, max_value):
        return max(min(num, max_value), min_value)

    def breakInfiniteLoop():
        set_value(name = "exit_angle", value = "Last computation failed.")
        delete_series(plot="nozzle_geometry",series="Upstream Curve")
        delete_series(plot="nozzle_geometry",series="Init. Exp. Curve")
        delete_series(plot="nozzle_geometry",series="Parabolic Curve")
        delete_series(plot="nozzle_geometry",series="Upstream Curve(M)")
        delete_series(plot="nozzle_geometry",series="Init. Exp. Curve(M)")
        delete_series(plot="nozzle_geometry",series="Parabolic Curve(M)")
        log_error("Infinite loop. Could not complete Run [" + str(calc_run_number) + "].", logger="Logs")
        return

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

    def exit_angle(Ln, iterative_radius, exit_radius):
        dy = exit_radius-iterative_radius
        delta_e = float(math.degrees(math.atan(dy/Ln)))
        return delta_e

    # First curve calculations

    theta_initial = math.degrees(math.asin(clamp((chamber_radius - (1.5 * throat_radius + throat_radius))/(1.5*throat_radius), -1, 1)))
    fc_steps = (-90-theta_initial)/20

    loop1_count = 0

    while True:
        x_fc = firstcurve_xfc(theta_FC = theta_FC, throat_radius = throat_radius)
        y_fc = firstcurve_yfc(theta_FC = theta_FC, throat_radius = throat_radius)
        first_curve_dict.update({ x_fc : y_fc })
        theta_FC = theta_FC + fc_steps
        loop1_count += 1

        if y_fc >= chamber_radius :
            theta_FC = theta_FC - 2*fc_steps
            break

        if loop1_count > 100000:
            breakInfiniteLoop()
            return

    # Second curve calculations

    # we have to use clamp() here because sometimes, due to floating point imprecision, we get a math domain error
    theta_SC = math.degrees(math.asin(clamp((throat_radius - (0.382 * throat_radius + throat_radius))/(0.382*throat_radius), -1, 1)))
    Sc_steps = (delta_n)/22

    loop2_count = 0
    while True:

        x_sc = secondcurve_xsc(theta_SC = theta_SC, throat_radius = throat_radius)
        y_sc = secondcurve_ysc(theta_SC = theta_SC, throat_radius = throat_radius)
        second_curve_dict.update({ x_sc : y_sc })
        theta_SC = theta_SC + Sc_steps
        loop2_count += 1

        if theta_SC >= (delta_n - 90) :
            theta_SC = theta_SC - 2*Sc_steps

            break

        if loop2_count > 100000:
            breakInfiniteLoop()
            return

    # Exit angle

    Ln = nozzle_length(throat_radius, expansion_ratio)
    delta_e = exit_angle(Ln, max(second_curve_dict.values()), exit_radius)

    # Update exit angle to display    
    set_value(name = "exit_angle", value = delta_e)
    # update value for export
    last_exit_angle = delta_e

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

    loop3_count = 0
    
    while True:
        px = a*second_curve_y**2+b*second_curve_y+c
        parabolic_x.append(px)
        parabolic_y.append(second_curve_y)
        second_curve_y = second_curve_y + parabolic_step
        loop3_count += 1
        
        if second_curve_y > exit_radius:
            break
        
        if loop3_count > 100000:
            breakInfiniteLoop()
            return

    ################################
    #   NOZZLE CALCS. END HERE
    ################################

    global last_first_curve, last_second_curve, last_parabolic_curve
    last_first_curve = [list(first_curve_dict.keys()), list(first_curve_dict.values())]
    last_second_curve = [list(second_curve_dict.keys()), list(second_curve_dict.values())]
    last_parabolic_curve = [parabolic_x, parabolic_y]

    #add data to plot
    add_line_series(name="Upstream Curve" , plot="nozzle_geometry",x=list(first_curve_dict.keys()), y=list(first_curve_dict.values()))
    add_line_series(name="Init. Exp. Curve" , plot="nozzle_geometry",x=list(second_curve_dict.keys()), y=list(second_curve_dict.values()))
    add_line_series(name="Parabolic Curve" , plot="nozzle_geometry",x=parabolic_x, y=parabolic_y)

    if get_value("mirror_check"):
        first_curve_y_mirrored = [-y for y in list(first_curve_dict.values())]
        second_curve_y_mirrored = [-y for y in list(second_curve_dict.values())]
        parabolic_y_mirrored = [-y for y in parabolic_y]
        
        add_line_series(name="Upstream Curve(M)" , plot="nozzle_geometry",x=list(first_curve_dict.keys()), y=first_curve_y_mirrored, color=(0,191,255,255))
        add_line_series(name="Init. Exp. Curve(M)" , plot="nozzle_geometry",x=list(second_curve_dict.keys()), y=second_curve_y_mirrored, color=(255,0,0,255))
        add_line_series(name="Parabolic Curve(M)" , plot="nozzle_geometry",x=parabolic_x, y=parabolic_y_mirrored, color=(127,255,0,255))
    else:
        #clean up mirrored curves (if any exists) so that they won't linger around after the checkbox is un-checked.
        delete_series(plot="nozzle_geometry",series="Upstream Curve(M)")
        delete_series(plot="nozzle_geometry",series="Init. Exp. Curve(M)")
        delete_series(plot="nozzle_geometry",series="Parabolic Curve(M)")

#FILE OPERATIONS BAR
with window("File I/O", width=960, height=60, no_close=True, no_move=True):
    set_window_pos("File I/O", 10, 10)
    add_input_text(name="filepath_field", label="Filepath", tip = "If the file is in the same directory with the script, you don't need\nto write the full path.")
    add_same_line()
    add_button("Import", callback = importFile)
    add_same_line()
    add_button("Export", callback = exportFile)

#INPUTS WINDOW
with window("Input", width=450, height=250, no_close=True):   
    set_window_pos("Input", 10, 80)
    add_text("Enter design parameters in float values.")
    add_spacing(count=6)
    add_input_text(name = "throat_radius_field", label = "Throat Radius (mm)")
    add_input_text(name = "exit_radius_field", label = "Exit Radius (mm)")
    add_input_text(name = "chamber_radius_field", label = "Chamber Radius (mm)")
    add_input_text(name = "delta_n_field", label = "Delta n (Throat Angle, deg)")
    add_input_text(name = "theta_FC_field", label = "Theta FC")
    add_spacing(count=6)
    add_button("Compute Nozzle Geometry", callback = computeNozzle)
    add_same_line()
    add_checkbox(name="mirror_check", label="Mirror Graph", tip="This only affects graph visualization. Exported data always uses the upper curve.")

#OUTPUTS WINDOW
with window("Output", width=500, height=560, no_close=True):
    set_window_pos("Output", 470, 80)
    add_input_text(name="exit_angle_output", label="Exit Angle (deg)", source="exit_angle", readonly=True, enabled=False)
    add_input_text(name="nozzle_expansion_ratio_output", label="Nozzle Expansion Ratio", source="nozzle_expansion_ratio", readonly=True, enabled=False)
    add_input_text(name="chamber_contraction_ratio_output", label="Chamber Contraction Ratio", source="chamber_contraction_ratio", readonly=True, enabled=False)
    add_plot(name="nozzle_geometry", label="Nozzle Geometry",
             x_axis_name="X Axis (mm)", y_axis_name = "Y Axis (mm)", equal_aspects = True)

#LOG WINDOW
with window("Log", width=450, height=300, no_close=True):
    set_window_pos("Log", 10, 340)
    add_logger("Logs", log_level=0, autosize_x = True, autosize_y = True)

start_dearpygui()
