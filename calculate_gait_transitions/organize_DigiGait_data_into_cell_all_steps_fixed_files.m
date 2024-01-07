% this function is used to ipload fixed stanc/swing traking sata for PVFlpO;Lbx1Cre;DTR 
% genotype (control and experimentals).

% output: swing_dur_all_steps_fixed_files, stance_dur_all_steps_fixed_files
% and cycle_dur_all_steps_fixed_files
% are sx2 cell arrays where cell location{1,1} and {1,2} contains values for control and
% experimentals info respectively. locations {2,1} and {2,2} contains mouse
% ID for each value in {1,1} and {1,2} respectively. 
% each cell is a 4x10,000x4 3D matrix where first dimention is the number of limbs
% (order: Left Forelimb, Right Forelimb, Left Hindlimb, Right Hindlimb). Second dimention is the number of step cycles (actual number
% is much lower, picked a large number to make matrix even (the rest of the
% elements are nans). third dimention is number of speeds (24,40,60,80). 
function organize_DigiGait_data_into_cell_all_steps_fixed_files

speed_vec = 20:20:80;
condition_types = ['C' , 'E'];
base_folder = 'V:\Behavior data\Data\DigiGait\PV project\PVFlpO;Lbx1Cre;DTR';
[~,~,raw] =xlsread([base_folder '\Copy of Paw area_all_animals_Mel2.0.xlsx']);
exclude = nan(length(raw),1);
speeds = nan(length(raw),1);
condition = strings(length(raw),1);
mouse_ID = nan(length(raw),1);

for i=2:length(raw)
    speeds(i) = raw{i,4};
    condition(i) = raw{i,6};
    exclude(i) = raw{i,10};
    fixed{i} = raw{i,12};
    mouse_ID(i) = raw{i,11};
end


swing_dur_all_steps_fixed_files = cell(2,2);
stance_dur_all_steps_fixed_files  = cell(2,2);
step_freq_all_steps_fixed_files  = cell(2,2);
stride_length_all_steps_fixed_files  = cell(2,2);
stance_to_swing_norm_all_steps_fixed_files  = cell(2,2);
swing_to_stance_norm_all_steps_fixed_files  = cell(2,2);


for i = 1:size(swing_dur_all_steps_fixed_files,1)
    for j = 1:size(swing_dur_all_steps_fixed_files,2) 
        swing_dur_all_steps_fixed_files{i,j} = nan(4,10000,4);
        stance_dur_all_steps_fixed_files{i,j} = nan(4,10000,4);
        step_freq_all_steps_fixed_files{i,j} = nan(4,10000,4);
        stride_length_all_steps_fixed_files{i,j} =  nan(4,10000,4);
        stance_to_swing_all_steps_fixed_files{i,j} = nan(4,10000,4);
        swing_to_stance_all_steps_fixed_files{i,j} = nan(4,10000,4);
    end
end

for i=1:length(condition_types)
    i
   for j = 1:length(speed_vec)
       j
        curr_speed = speed_vec(j);
        curr_ind = find(condition==condition_types(i) & speeds==curr_speed & exclude==0 & strcmp(fixed, 'fixed')');    
            for m=1:length(curr_ind)                          
                m
                curr_dob= raw{curr_ind(m),3};
                slash_ind = find(curr_dob=='/');
                animal_name = raw{curr_ind(m),1};
                stance_file_name = [base_folder '\PVFlpO;Lbx1Cre;DTR;Ai65 ',curr_dob(1:slash_ind(1)-1),'-',curr_dob(slash_ind(1)+1:slash_ind(2)-1),'-20\Paw contact area new\stance_inds_',animal_name,'_',num2str(curr_speed),'_cms_Mel'];
                swing_file_name = [base_folder '\PVFlpO;Lbx1Cre;DTR;Ai65 ',curr_dob(1:slash_ind(1)-1),'-',curr_dob(slash_ind(1)+1:slash_ind(2)-1),'-20\Paw contact area new\swing_inds_',animal_name,'_',num2str(curr_speed),'_cms_Mel'];
                load(stance_file_name)
                load(swing_file_name)   
                [swing_dur{m}, stance_dur{m}, step_freq{m}, stance_to_swing_norm{m}, swing_to_stance_norm{m}] = get_locmotion_param_for_fixed_files(stance_ind_start_corr_corr, swing_ind_start_corr_corr);
                curr_mouse_ID = mouse_ID(curr_ind(m));
                for k = 1:4
                    curr_mat_swing =  swing_dur{m}{k}';
                    last_ind_swing =  find(isnan(swing_dur_all_steps_fixed_files{1,i}(k,:,j)),1);
                    swing_dur_all_steps_fixed_files{1,i}(k,last_ind_swing:last_ind_swing+length(curr_mat_swing)-1,j) = curr_mat_swing;
                    swing_dur_all_steps_fixed_files{2,i}(k,last_ind_swing:last_ind_swing+length(curr_mat_swing)-1,j) =  curr_mouse_ID ;

                    curr_mat_stance =  stance_dur{m}{k}';
                    last_ind_stance =  find(isnan(stance_dur_all_steps_fixed_files{1,i}(k,:,j)),1);
                    stance_dur_all_steps_fixed_files{1,i}(k,last_ind_stance:last_ind_stance+length(curr_mat_stance)-1,j) = curr_mat_stance;
                    stance_dur_all_steps_fixed_files{2,i}(k,last_ind_stance:last_ind_stance+length(curr_mat_stance)-1,j) =  curr_mouse_ID ;

                    curr_mat_freq =  step_freq{m}{k}';
                    last_ind_freq =  find(isnan(step_freq_all_steps_fixed_files{1,i}(k,:,j)),1);
                    step_freq_all_steps_fixed_files{1,i}(k,last_ind_freq:last_ind_freq+length(curr_mat_freq)-1,j) = curr_mat_freq; 
                    step_freq_all_steps_fixed_files{2,i}(k,last_ind_freq:last_ind_freq+length(curr_mat_freq)-1,j) =  curr_mouse_ID ;

                    curr_mat_stride = (curr_mat_stance+curr_mat_swing)*curr_speed; %calculate stride length in cm. 
                    last_ind_stride_length =  find(isnan(step_freq_all_steps_fixed_files{1,i}(k,:,j)),1);
                    stride_length_all_steps_fixed_files{1,i}(k,last_ind_stride_length:last_ind_stride_length+length(curr_mat_stride)-1,j) = curr_mat_stride; 
                    stride_length_all_steps_fixed_files{2,i}(k,last_ind_stride_length:last_ind_stride_length+length(curr_mat_stride)-1,j) =  curr_mouse_ID ;

                    curr_mat_stance_to_swing =  stance_to_swing_norm{m}{k}';
                    last_ind_stance_to_swing =  find(isnan(stance_to_swing_all_steps_fixed_files{1,i}(k,:,j)),1);
                    stance_to_swing_all_steps_fixed_files{1,i}(k,last_ind_stance_to_swing:last_ind_stance_to_swing+length(curr_mat_stance_to_swing)-1,j) = curr_mat_stance_to_swing; 
                    stance_to_swing_all_steps_fixed_files{2,i}(k,last_ind_stance_to_swing:last_ind_stance_to_swing+length(curr_mat_stance_to_swing)-1,j) =  curr_mouse_ID ;

                    curr_mat_swing_to_stance =  swing_to_stance_norm{m}{k}';
                    last_ind_swing_to_stance =  find(isnan(swing_to_stance_all_steps_fixed_files{1,i}(k,:,j)),1);
                    swing_to_stance_all_steps_fixed_files{1,i}(k,last_ind_swing_to_stance:last_ind_swing_to_stance+length(curr_mat_swing_to_stance)-1,j) = curr_mat_swing_to_stance; 
                    swing_to_stance_all_steps_fixed_files{2,i}(k,last_ind_swing_to_stance:last_ind_swing_to_stance+length(curr_mat_swing_to_stance)-1,j) =  curr_mouse_ID ;

                end
                clear curr_mat* *corr_corr
            end
            clear swing_dur stance_dur step_freq_dur
   end
end



% save('V:\Nofar\Behavior\DigiGait\data in m files\swing_dur_all_steps_fixed_files','swing_dur_all_steps_fixed_files')
% save('V:\Nofar\Behavior\DigiGait\data in m files\stance_dur_all_steps_fixed_files','stance_dur_all_steps_fixed_files')
% save('V:\Nofar\Behavior\DigiGait\data in m files\step_freq_all_steps_fixed_files','step_freq_all_steps_fixed_files')
% save('V:\Nofar\Behavior\DigiGait\data in m files\stance_to_swing_all_steps_fixed_files','stance_to_swing_all_steps_fixed_files')
% save('V:\Nofar\Behavior\DigiGait\data in m files\swing_to_stance_all_steps_fixed_files','swing_to_stance_all_steps_fixed_files')

end
