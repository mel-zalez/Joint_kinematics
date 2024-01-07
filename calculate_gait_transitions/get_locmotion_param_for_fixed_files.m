function [swing_dur, stance_dur, step_freq, stance_to_swing_norm, swing_to_stance_norm] = get_locmotion_param_for_fixed_files(stance_ind_start_corr_corr, swing_ind_start_corr_corr)

frame_rate = 165;


for i= 1:4
    num_step_cycle = min([length(stance_ind_start_corr_corr{1,i})-1 length(swing_ind_start_corr_corr{1,i})-1]);
    for j = 1: num_step_cycle
        curr_stance_st = stance_ind_start_corr_corr{1,i}(j);
        next_stance_st = stance_ind_start_corr_corr{1,i}(j+1);
        curr_swing_st = swing_ind_start_corr_corr{1,i}(swing_ind_start_corr_corr{1,i}> curr_stance_st &  swing_ind_start_corr_corr{1,i}< next_stance_st);
        
        stance_dur{i}(j) = (curr_swing_st-curr_stance_st)/frame_rate;
        curr_step_cycle_dur = (next_stance_st-curr_stance_st)/frame_rate;
        swing_dur{i}(j) = curr_step_cycle_dur-stance_dur{i}(j);
        stance_to_swing_norm{i}(j) = (curr_swing_st - curr_stance_st) / (next_stance_st - curr_stance_st);% normalized stance to swing
        step_freq{i}(j) = 1./curr_step_cycle_dur;
    end
end
    
%     to calculate normalized swing to stance transition, define step cycle
%     from begining of one swing to next

for i= 1:4
    num_step_cycle = min([length(stance_ind_start_corr_corr{1,i})-1 length(swing_ind_start_corr_corr{1,i})-1]);
    for j = 1: num_step_cycle
        curr_swing_st = swing_ind_start_corr_corr{1,i}(j);
        next_swing_st = swing_ind_start_corr_corr{1,i}(j+1);
        curr_stance_st = stance_ind_start_corr_corr{1,i}(stance_ind_start_corr_corr{1,i}> curr_swing_st &  stance_ind_start_corr_corr{1,i}< next_swing_st);
        swing_to_stance_norm{i}(j) = (curr_stance_st - curr_swing_st) / ( next_swing_st - curr_swing_st);% normalized stance to swing
    end
end