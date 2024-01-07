% each cell is a 4x10,000x4 3D matrix where first dimension is the number of limbs
% (order: Left Forelimb, Right Forelimb, Left Hindlimb, Right Hindlimb). Second dimention is the number of step cycles (actual number
% is much lower, picked a large number to make matrix even (the rest of the
% elements are nans). third dimension is number of speeds (20, 40, 60, 80). 


function [p_swing_corr, p_stance_corr,p_freq_corr, p_stride, p_st_to_sw_corr, p_sw_to_st_corr]  = compare_locomotion_parameters_all_steps_fixed_files()

load('V:\Nofar\Behavior\DigiGait\data in m files\swing_dur_all_steps_fixed_files.mat')
load('V:\Nofar\Behavior\DigiGait\data in m files\stance_dur_all_steps_fixed_files.mat')
load('V:\Nofar\Behavior\DigiGait\data in m files\step_freq_all_steps_fixed_files.mat')
load('V:\Nofar\Behavior\DigiGait\data in m files\stride_length_all_steps_fixed_files.mat')
load('V:\Nofar\Behavior\DigiGait\data in m files\stance_to_swing_all_steps_fixed_files.mat')
load('V:\Nofar\Behavior\DigiGait\data in m files\swing_to_stance_all_steps_fixed_files.mat')
speed_vec = 20:20:80; 


for i = 1:length(speed_vec)
    curr_ctr_swing_mat = [swing_dur_all_steps_fixed_files{1,1}(3,:,i) swing_dur_all_steps_fixed_files{1,1}(4,:,i)];
    curr_ctr_swing_mat = curr_ctr_swing_mat(~isnan(curr_ctr_swing_mat));
%     length(  curr_ctr_swing_mat)
    mean_swing_ctr(i) = mean(curr_ctr_swing_mat);
    sem_swing_ctr(i) = std(curr_ctr_swing_mat)/sqrt(length(curr_ctr_swing_mat));
    curr_ctr_swing_ID = [swing_dur_all_steps_fixed_files{2,1}(3,:,i) swing_dur_all_steps_fixed_files{2,1}(4,:,i)];
    curr_ctr_swing_ID = curr_ctr_swing_ID(~isnan(curr_ctr_swing_ID));
    
    curr_exp_swing_mat = [swing_dur_all_steps_fixed_files{1,2}(3,:,i) swing_dur_all_steps_fixed_files{1,2}(4,:,i)]; 
    curr_exp_swing_mat = curr_exp_swing_mat(~isnan(curr_exp_swing_mat));
    mean_swing_exp(i) = mean(curr_exp_swing_mat);
    sem_swing_exp(i) = std(curr_exp_swing_mat)/sqrt(length(curr_exp_swing_mat));
    curr_exp_swing_ID = [swing_dur_all_steps_fixed_files{2,2}(3,:,i) swing_dur_all_steps_fixed_files{2,2}(4,:,i)];
    curr_exp_swing_ID =  curr_exp_swing_ID(~isnan(curr_exp_swing_ID));
       
    swing_vals = [curr_ctr_swing_mat';curr_exp_swing_mat'];
    swing_dummy_vars = [zeros(length(curr_ctr_swing_mat),1);1+zeros(length(curr_exp_swing_mat),1)];
    swing_mouse_id = [curr_ctr_swing_ID';curr_exp_swing_ID'];
    T_swing = table( swing_vals,swing_dummy_vars,swing_mouse_id,'VariableNames',{'Swing_time','MouseType','MouseID'});
%     lme_swing = fitlme( T_swing,'Swing_time~1+MouseType+(1|MouseID)+(MouseType|MouseID)');
%      lme_swing = fitlme( T_swing,'Swing_time~1+MouseType +(MouseType|MouseID)');
    lme_swing = fitlme( T_swing,'Swing_time~1+MouseType+(1|MouseID)');
    p_swing(i) = lme_swing.Coefficients.pValue(2);
    
%     
    curr_ctr_stance_mat = [stance_dur_all_steps_fixed_files{1,1}(3,:,i) stance_dur_all_steps_fixed_files{1,1}(4,:,i)];
    curr_ctr_stance_mat = curr_ctr_stance_mat(~isnan(curr_ctr_stance_mat));
    mean_stance_ctr(i) = mean(curr_ctr_stance_mat);
    sem_stance_ctr(i) = std(curr_ctr_stance_mat)/sqrt(length(curr_ctr_stance_mat));
    curr_ctr_stance_ID = [stance_dur_all_steps_fixed_files{2,1}(3,:,i) stance_dur_all_steps_fixed_files{2,1}(4,:,i)];
    curr_ctr_stance_ID = curr_ctr_stance_ID(~isnan(curr_ctr_stance_ID));
    
    curr_exp_stance_mat = [stance_dur_all_steps_fixed_files{1,2}(3,:,i) stance_dur_all_steps_fixed_files{1,2}(4,:,i)]; 
    curr_exp_stance_mat = curr_exp_stance_mat(~isnan(curr_exp_stance_mat));
    mean_stance_exp(i) = mean(curr_exp_stance_mat);
    sem_stance_exp(i) = std(curr_exp_stance_mat)/sqrt(length(curr_exp_stance_mat));
    curr_exp_stance_ID = [stance_dur_all_steps_fixed_files{2,2}(3,:,i) stance_dur_all_steps_fixed_files{2,2}(4,:,i)];
    curr_exp_stance_ID =  curr_exp_stance_ID(~isnan(curr_exp_stance_ID));
       
    stance_vals = [curr_ctr_stance_mat';curr_exp_stance_mat'];
    stance_dummy_vars = [zeros(length(curr_ctr_stance_mat),1);1+zeros(length(curr_exp_stance_mat),1)];
    stance_mouse_id = [curr_ctr_stance_ID';curr_exp_stance_ID'];
    T_stance = table(stance_vals,stance_dummy_vars,stance_mouse_id,'VariableNames',{'Stance_time','MouseType','MouseID'});
%     lme_stance = fitlme(T_stance,'Stance_time~1+MouseType+(1|MouseID)+(MouseType|MouseID)');
%     lme_stance = fitlme(T_stance,'Stance_time~1+MouseType+(MouseType|MouseID)');
    lme_stance = fitlme(T_stance,'Stance_time~1+MouseType+(1|MouseID)');
    p_stance(i) = lme_stance.Coefficients.pValue(2);
    
%     


% stride duration:    
   
    curr_ctr_stride_time_mat = [stance_dur_all_steps_fixed_files{1,1}(3,:,i) stance_dur_all_steps_fixed_files{1,1}(4,:,i)];
    mean_stance_ctr(i) = mean(curr_ctr_stance_mat);
    sem_stance_ctr(i) = std(curr_ctr_stance_mat)/sqrt(length(curr_ctr_stance_mat));
    curr_ctr_stance_ID = [stance_dur_all_steps_fixed_files{2,1}(3,:,i) stance_dur_all_steps_fixed_files{2,1}(4,:,i)];
    curr_ctr_stance_ID = curr_ctr_stance_ID(~isnan(curr_ctr_stance_ID));
    
    curr_exp_stance_mat = [stance_dur_all_steps_fixed_files{1,2}(3,:,i) stance_dur_all_steps_fixed_files{1,2}(4,:,i)]; 
    curr_exp_stance_mat = curr_exp_stance_mat(~isnan(curr_exp_stance_mat));
    mean_stance_exp(i) = mean(curr_exp_stance_mat);
    sem_stance_exp(i) = std(curr_exp_stance_mat)/sqrt(length(curr_exp_stance_mat));
    curr_exp_stance_ID = [stance_dur_all_steps_fixed_files{2,2}(3,:,i) stance_dur_all_steps_fixed_files{2,2}(4,:,i)];
    curr_exp_stance_ID =  curr_exp_stance_ID(~isnan(curr_exp_stance_ID));
       
    stance_vals = [curr_ctr_stance_mat';curr_exp_stance_mat'];
    stance_dummy_vars = [zeros(length(curr_ctr_stance_mat),1);1+zeros(length(curr_exp_stance_mat),1)];
    stance_mouse_id = [curr_ctr_stance_ID';curr_exp_stance_ID'];
    T_stance = table(stance_vals,stance_dummy_vars,stance_mouse_id,'VariableNames',{'Stance_time','MouseType','MouseID'});
%     lme_stance = fitlme(T_stance,'Stance_time~1+MouseType+(1|MouseID)+(MouseType|MouseID)');
%     lme_stance = fitlme(T_stance,'Stance_time~1+MouseType+(MouseType|MouseID)');
    lme_stance = fitlme(T_stance,'Stance_time~1+MouseType+(1|MouseID)');
    p_stance(i) = lme_stance.Coefficients.pValue(2);
    
%     






















    curr_ctr_freq_mat = [step_freq_all_steps_fixed_files{1,1}(3,:,i) step_freq_all_steps_fixed_files{1,1}(4,:,i)];
    curr_ctr_freq_mat = curr_ctr_freq_mat(~isnan(curr_ctr_freq_mat));
    mean_freq_ctr(i) = mean(curr_ctr_freq_mat);
    sem_freq_ctr(i) = std(curr_ctr_freq_mat)/sqrt(length(curr_ctr_freq_mat));
    curr_ctr_freq_ID = [step_freq_all_steps_fixed_files{2,1}(3,:,i) step_freq_all_steps_fixed_files{2,1}(4,:,i)];
    curr_ctr_freq_ID = curr_ctr_freq_ID(~isnan(curr_ctr_freq_ID));
    
    curr_exp_freq_mat = [step_freq_all_steps_fixed_files{1,2}(3,:,i) step_freq_all_steps_fixed_files{1,2}(4,:,i)]; 
    curr_exp_freq_mat = curr_exp_freq_mat(~isnan(curr_exp_freq_mat));
    mean_freq_exp(i) = mean(curr_exp_freq_mat);
    sem_freq_exp(i) = std(curr_exp_freq_mat)/sqrt(length(curr_exp_freq_mat));
    curr_exp_freq_ID = [step_freq_all_steps_fixed_files{2,2}(3,:,i) step_freq_all_steps_fixed_files{2,2}(4,:,i)];
    curr_exp_freq_ID =  curr_exp_freq_ID(~isnan(curr_exp_freq_ID));
       
    freq_vals = [curr_ctr_freq_mat';curr_exp_freq_mat'];
    freq_dummy_vars = [zeros(length(curr_ctr_freq_mat),1);1+zeros(length(curr_exp_freq_mat),1)];
    freq_mouse_id = [curr_ctr_freq_ID';curr_exp_freq_ID'];
    T_freq = table(freq_vals,freq_dummy_vars,freq_mouse_id,'VariableNames',{'Step_freq','MouseType','MouseID'});
%     lme_freq = fitlme(T_freq,'Step_freq~1+MouseType+(1|MouseID)+(MouseType|MouseID)');
%     lme_freq = fitlme(T_freq,'Step_freq~1+MouseType+(MouseType|MouseID)');
    lme_freq = fitlme(T_freq,'Step_freq~1+MouseType+(1|MouseID)');
    p_freq(i) = lme_freq.Coefficients.pValue(2);
    
    
%     
    curr_ctr_stride_mat = [stride_length_all_steps_fixed_files{1,1}(3,:,i) stride_length_all_steps_fixed_files{1,1}(4,:,i)];
    curr_ctr_stride_mat =  curr_ctr_stride_mat(~isnan(curr_ctr_stride_mat));
    mean_stride_ctr(i) = mean(curr_ctr_stride_mat);
    sem_stride_ctr(i) = std(curr_ctr_stride_mat)/sqrt(length(curr_ctr_stride_mat));
    curr_ctr_stride_ID = [stride_length_all_steps_fixed_files{2,1}(3,:,i) stride_length_all_steps_fixed_files{2,1}(4,:,i)];
    curr_ctr_stride_ID = curr_ctr_stride_ID(~isnan(curr_ctr_stride_ID));
    
    curr_exp_stride_mat = [stride_length_all_steps_fixed_files{1,2}(3,:,i) stride_length_all_steps_fixed_files{1,2}(4,:,i)]; 
    curr_exp_stride_mat = curr_exp_stride_mat(~isnan(curr_exp_stride_mat));
    mean_stride_exp(i) = mean(curr_exp_stride_mat);
    sem_stride_exp(i) = std(curr_exp_stride_mat)/sqrt(length(curr_exp_stride_mat));
    curr_exp_stride_ID = [stride_length_all_steps_fixed_files{2,2}(3,:,i) stride_length_all_steps_fixed_files{2,2}(4,:,i)];
    curr_exp_stride_ID =  curr_exp_stride_ID(~isnan(curr_exp_stride_ID));
       
    stride_vals = [curr_ctr_stride_mat';curr_exp_stride_mat'];
    stride_dummy_vars = [zeros(length(curr_ctr_stride_mat),1);1+zeros(length(curr_exp_stride_mat),1)];
    stride_mouse_id = [curr_ctr_stride_ID';curr_exp_stride_ID'];
    T_stride = table(stride_vals,stride_dummy_vars,stride_mouse_id,'VariableNames',{'Stride_length','MouseType','MouseID'});
%     lme_freq = fitlme(T_freq,'Step_freq~1+MouseType+(1|MouseID)+(MouseType|MouseID)');
%     lme_freq = fitlme(T_freq,'Step_freq~1+MouseType+(MouseType|MouseID)');
    lme_stride_length = fitlme(T_stride,'Stride_length~1+MouseType+(1|MouseID)');
    p_stride(i) = lme_stride_length.Coefficients.pValue(2);
    
%     
    curr_ctr_st_to_sw_mat = [stance_to_swing_all_steps_fixed_files{1,1}(3,:,i) stance_to_swing_all_steps_fixed_files{1,1}(4,:,i)];
    curr_ctr_st_to_sw_mat =  curr_ctr_st_to_sw_mat(~isnan(curr_ctr_st_to_sw_mat));
    mean_st_to_sw_ctr(i) = mean(curr_ctr_st_to_sw_mat);
    sem_st_to_sw_ctr(i) = std(curr_ctr_st_to_sw_mat)/sqrt(length(curr_ctr_st_to_sw_mat));
    curr_ctr_st_to_sw_ID = [stance_to_swing_all_steps_fixed_files{2,1}(3,:,i) stance_to_swing_all_steps_fixed_files{2,1}(4,:,i)];
    curr_ctr_st_to_sw_ID = curr_ctr_st_to_sw_ID(~isnan(curr_ctr_st_to_sw_ID));
    
    curr_exp_st_to_sw_mat = [stance_to_swing_all_steps_fixed_files{1,2}(3,:,i) stance_to_swing_all_steps_fixed_files{1,2}(4,:,i)]; 
    curr_exp_st_to_sw_mat =  curr_exp_st_to_sw_mat(~isnan(curr_exp_st_to_sw_mat));
    mean_st_to_sw_exp(i) = mean(curr_exp_st_to_sw_mat);
    sem_st_to_sw_exp(i) = std(curr_exp_st_to_sw_mat)/sqrt(length(curr_exp_st_to_sw_mat));
    curr_exp_st_to_sw_ID = [stance_to_swing_all_steps_fixed_files{2,2}(3,:,i) stance_to_swing_all_steps_fixed_files{2,2}(4,:,i)];
    curr_exp_st_to_sw_ID =  curr_exp_st_to_sw_ID(~isnan(curr_exp_st_to_sw_ID));
       
    st_to_sw_vals = [curr_ctr_st_to_sw_mat';curr_exp_st_to_sw_mat'];
    st_to_sw_dummy_vars = [zeros(length(curr_ctr_st_to_sw_mat),1);1+zeros(length(curr_exp_st_to_sw_mat),1)];
    st_to_sw_mouse_id = [curr_ctr_st_to_sw_ID';curr_exp_st_to_sw_ID'];
    T_st_to_sw = table(st_to_sw_vals,st_to_sw_dummy_vars,st_to_sw_mouse_id,'VariableNames',{'st_to_sw','MouseType','MouseID'});
%     lme_freq = fitlme(T_freq,'Step_freq~1+MouseType+(1|MouseID)+(MouseType|MouseID)');
%     lme_freq = fitlme(T_freq,'Step_freq~1+MouseType+(MouseType|MouseID)');
    lme_st_to_sw = fitlme(T_st_to_sw,'st_to_sw~1+MouseType+(1|MouseID)');
    p_st_to_sw(i) = lme_st_to_sw.Coefficients.pValue(2);
    
 %     
    curr_ctr_sw_to_st_mat = [swing_to_stance_all_steps_fixed_files{1,1}(3,:,i) swing_to_stance_all_steps_fixed_files{1,1}(4,:,i)];
    curr_ctr_sw_to_st_mat =  curr_ctr_sw_to_st_mat(~isnan(curr_ctr_sw_to_st_mat));
    mean_sw_to_st_ctr(i) = mean(curr_ctr_sw_to_st_mat);
    sem_sw_to_st_ctr(i) = std(curr_ctr_sw_to_st_mat)/sqrt(length(curr_ctr_sw_to_st_mat));
    curr_ctr_sw_to_st_ID = [swing_to_stance_all_steps_fixed_files{2,1}(3,:,i) swing_to_stance_all_steps_fixed_files{2,1}(4,:,i)];

    curr_ctr_sw_to_st_ID = curr_ctr_sw_to_st_ID(~isnan(curr_ctr_sw_to_st_ID));
    
    curr_exp_sw_to_st_mat = [swing_to_stance_all_steps_fixed_files{1,2}(3,:,i) swing_to_stance_all_steps_fixed_files{1,2}(4,:,i)]; 
    curr_exp_sw_to_st_mat =  curr_exp_sw_to_st_mat(~isnan(curr_exp_sw_to_st_mat));
    mean_sw_to_st_exp(i) = mean(curr_exp_sw_to_st_mat);
    sem_sw_to_st_exp(i) = std(curr_exp_sw_to_st_mat)/sqrt(length(curr_exp_sw_to_st_mat));
    curr_exp_sw_to_st_ID = [swing_to_stance_all_steps_fixed_files{2,2}(3,:,i) swing_to_stance_all_steps_fixed_files{2,2}(4,:,i)];

    curr_exp_sw_to_st_ID =  curr_exp_sw_to_st_ID(~isnan(curr_exp_sw_to_st_ID));
       
    sw_to_st_vals = [curr_ctr_sw_to_st_mat';curr_exp_sw_to_st_mat'];
    sw_to_st_dummy_vars = [zeros(length(curr_ctr_sw_to_st_mat),1);1+zeros(length(curr_exp_sw_to_st_mat),1)];
    sw_to_st_mouse_id = [curr_ctr_sw_to_st_ID';curr_exp_sw_to_st_ID'];
    T_sw_to_st = table(sw_to_st_vals,sw_to_st_dummy_vars,sw_to_st_mouse_id,'VariableNames',{'sw_to_st','MouseType','MouseID'});
%     lme_freq = fitlme(T_freq,'Step_freq~1+MouseType+(1|MouseID)+(MouseType|MouseID)');
%     lme_freq = fitlme(T_freq,'Step_freq~1+MouseType+(MouseType|MouseID)');
    lme_sw_to_st = fitlme(T_sw_to_st,'sw_to_st~1+MouseType+(1|MouseID)');
    p_sw_to_st(i) = lme_sw_to_st.Coefficients.pValue(2);
       
    
%     [~,p_corr(1:3,i)] = find_holmbonferroni([p_swing(i),p_stance(i),p_freq(i)],0.05);%correct for multiple comparisons 
end

[~,p_swing_corr] = find_holmbonferroni(p_swing,0.05);%correct for multiple comparisons
[~,p_freq_corr] = find_holmbonferroni(p_freq,0.05);
[~,p_stance_corr] = find_holmbonferroni(p_stance,0.05);
[~,p_stride_corr] = find_holmbonferroni(p_stride,0.05);
[~,p_st_to_sw_corr] = find_holmbonferroni(p_st_to_sw,0.05);
[~,p_sw_to_st_corr] = find_holmbonferroni(p_sw_to_st,0.05);



figure
hold on
errorbar(mean_swing_ctr,sem_swing_ctr,'k');
errorbar(mean_swing_exp,sem_swing_exp,'r');

  
for i = 1:length(speed_vec)
    if p_swing_corr(i)<0.05 
        max_mean = max([mean_swing_ctr(i) mean_swing_exp(i)]);
        max_sem = max([sem_swing_ctr(i) sem_swing_exp(i)]);
        aster_loc = max_mean+max_sem+0.005;
        plot(i,aster_loc,'*','color','k')
    end
end
title('swing time')
% xticklabels({'5';'10';'20';'30';'40';'50';'60';'70';'80'})
xticks(1:1:4);
xticklabels({'20';'40';'60';'80'})
xlim([0.6 4.4])

set(gca, 'TickDir', 'out')
% saveas(gcf,'V:\Nofar\Behavior\DigiGait\Analysis\swing_time_hindlimbs_all_steps.emf')

figure
hold on
errorbar(mean_stance_ctr,sem_stance_ctr,'k');
errorbar(mean_stance_exp,sem_stance_exp,'r');
for i = 1:length(speed_vec)
    if p_stance_corr(i)<0.05 
        max_mean = max([mean_stance_ctr(i) mean_stance_exp(i)]);
        max_sem = max([sem_stance_ctr(i) sem_stance_exp(i)]);
        aster_loc = max_mean+max_sem+0.005;
        plot(i,aster_loc,'*','color','k')
    end
end
title('stance time')
% xticklabels({'5';'10';'20';'30';'40';'50';'60';'70';'80'})
set(gca, 'TickDir', 'out')
xticks(1:1:4);
xticklabels({'20';'40';'60';'80'})
xlim([0.6 4.4])
ylim([0.035 0.16]);shg
% saveas(gcf,'V:\Nofar\Behavior\DigiGait\Analysis\stance_time_hindlimbs_all_steps.emf')


figure
hold on
errorbar(mean_freq_ctr,sem_freq_ctr,'k');
errorbar(mean_freq_exp,sem_freq_exp,'r');
for i = 1:length(speed_vec)
    if p_freq_corr(i)<0.05 
        max_mean = max([mean_freq_ctr(i) mean_freq_exp(i)]);
        max_sem = max([sem_freq_ctr(i) sem_freq_exp(i)]);
        aster_loc = max_mean+max_sem+0.005;
        plot(i,aster_loc,'*','color','k')
    end
end
title('step frequency')
% xticklabels({'5';'10';'20';'30';'40';'50';'60';'70';'80'})
xticks(1:1:4);
xticklabels({'20';'40';'60';'80'})
xlim([0.6 4.4])
set(gca, 'TickDir', 'out')
% saveas(gcf,'V:\Nofar\Behavior\DigiGait\Analysis\step_frequency_hindlimbs_all_steps.emf')

figure
hold on
errorbar(mean_stride_ctr,sem_stride_ctr,'k');
errorbar(mean_stride_exp,sem_stride_exp,'r');
for i = 1:length(speed_vec)
    if p_stride_corr(i)<0.05 
        max_mean = max([mean_stride_ctr(i) mean_stride_exp(i)]);
        max_sem = max([sem_stride_ctr(i) sem_stride_exp(i)]);
        aster_loc = max_mean+max_sem+0.005;
        plot(i,aster_loc,'*','color','k')
    end
end
title('stride length')
% xticklabels({'5';'10';'20';'30';'40';'50';'60';'70';'80'})
xticks(1:1:4);
xticklabels({'20';'40';'60';'80'})
xlim([0.6 4.4])
ylim([4 9]);shg
set(gca, 'TickDir', 'out')
% saveas(gcf,'V:\Nofar\Behavior\DigiGait\Analysis\stride_length_hindlimbs_all_steps.emf')

figure
hold on
errorbar(mean_st_to_sw_ctr,sem_st_to_sw_ctr,'k');
errorbar(mean_st_to_sw_exp,sem_st_to_sw_exp,'r');

for i = 1:length(speed_vec)
    if p_st_to_sw_corr(i)<0.05 
        max_mean = max([mean_st_to_sw_ctr(i) mean_st_to_sw_exp(i)]);
        max_sem = max([sem_st_to_sw_ctr(i) sem_st_to_sw_exp(i)]);
        aster_loc = max_mean+max_sem+0.005;
        plot(i,aster_loc,'*','color','k')
    end
end
title('stance to swing transition phase')
% xticklabels({'5';'10';'20';'30';'40';'50';'60';'70';'80'})
xticks(1:1:4);
xticklabels({'20';'40';'60';'80'})
xlim([0.6 4.4])
ylim([0.35 0.75]);shg
set(gca, 'TickDir', 'out')
% saveas(gcf,'V:\Nofar\Behavior\DigiGait\Analysis\stance_to_swing_transition_hindlimbs_all_steps.emf')

figure
hold on
errorbar(mean_sw_to_st_ctr,sem_sw_to_st_ctr,'k');
errorbar(mean_sw_to_st_exp,sem_sw_to_st_exp,'r');

for i = 1:length(speed_vec)
    if p_sw_to_st_corr(i)<0.05 
        max_mean = max([mean_sw_to_st_ctr(i) mean_sw_to_st_exp(i)]);
        max_sem = max([sem_sw_to_st_ctr(i) sem_sw_to_st_exp(i)]);
        aster_loc = max_mean+max_sem+0.005;
        plot(i,aster_loc,'*','color','k')
    end
end
title('swing to stance transition phase')
% xticklabels({'5';'10';'20';'30';'40';'50';'60';'70';'80'})
xticks(1:1:4);
xticklabels({'20';'40';'60';'80'})
xlim([0.6 4.4])
ylim([0.25 0.65]);shg
set(gca, 'TickDir', 'out')
% saveas(gcf,'V:\Nofar\Behavior\DigiGait\Analysis\swing_to_stance_transition_hindlimbs_all_steps.emf')

end

