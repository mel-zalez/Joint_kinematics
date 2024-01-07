# Joint_kinematics

# calculate joint angles

Following DeepLabCut joint tracking, fill out video information in Excel file and use generate_new_skeleton.m to triangulate the position of the knee given the coordinates of the hip and ankle, the lengths of the tibia and femur in millimeteres, and the DeepLabCut tracking of the pixels to millimeters conversion grid. Use find_swing_stance_time.m to phase flag the stance and swing indices of each step cycle. Use plot_angle_vs_step_cycle_DTR.m to calculate joint angles for individual mice and compare_angles_DTR_all_steps to calculate joint angles and compare between the experimental and control mice. Use compare_locomotion_parameters_all_steps_files.m to calculate stride frequency, length, and duration between the experimental and control mice at varied speeds.

# calculate gait transitions

Following forelimb and hindlimb paw tracking with Digigait software, use get_locomotion_param_for_filed_files.m to extract stance and swing indices per paw per step cycle. Fill out video information in Excel file and use organize_DigiGait_data_into_cell_all_steps_fixed_files.m to organize stance and swing indeces for experimental and control mice at varied speeds and use plot_gait_transition_new_fixed_files.m to plot gait transitions with respective attraction and persistance. 
