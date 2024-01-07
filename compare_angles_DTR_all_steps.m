% ctr = {'EN0_2.17.21_Bl_M_20cms'; 'EN0_2.17.21_Wh_M_20cms'; 'EN1_2.17.21_Wh_F_20cms'};
% exp = {'EN0_2.16.21_Br_M_20cms'; 'EN1_2.16.21_Br_M_20cms';'EN0_2.17.21_Br_F';'EN0_2.17.21_Wh_F_20cms'};
function [p_AEP, p_PEP] = compare_angles_DTR_all_steps(ctr, exp)
exp_name = 'PV_Project_DTR';

% call plot_angle_vs_step_cycle_DTR for controls and experimentals to get: 
% mean/std_angle_ctr/exp is a nx4x1001 matix where n is the number of animals,
% 4 is the number of angles(hip, knee, ankle, mtp) and 1001 is the
% time axis). 
%  mean/std_swing_ctr/exp is 1xn matrice where n is the number of animals. 
% mean/std_AEP/PEP_ctr/exp is a 4xn matrix where 4 is the number of angles and 
% n is the # of animals 

for i = 1:size(ctr,1)
    [angle_by_cycle_ctr{i}, swing_ctr{i}, AEP_ctr{i},PEP_ctr{i},cycle_dur_ctr{i}, stance_dur_ctr{i}, swing_dur_ctr{i},step_num_ctr(i)] =  plot_angle_vs_step_cycle_DTR(exp_name,ctr{i});
    curr_num_trials = size(angle_by_cycle_ctr{i},2);
    ctr_ID{i}= i*ones(curr_num_trials,1);
end
all_angles_by_cycle_ctr = cat(2,angle_by_cycle_ctr{:});
ctr_ID = cat(1,ctr_ID{:});
all_swing_ctr = cat(2,swing_ctr{:});
all_AEP_ctr = cat(2,AEP_ctr{:});
all_PEP_ctr = cat(2,PEP_ctr{:});
N_ctr = size(all_angles_by_cycle_ctr,2);


for i = 1:size(exp,1)
    [angles_by_cycle_exp{i}, swing_exp{i},AEP_exp{i},PEP_exp{i},cycle_dur_exp{i}, stance_dur_exp{i}, swing_dur_exp{i},step_num_exp(i)] =  plot_angle_vs_step_cycle_DTR(exp_name,exp{i});
    curr_num_trials = size(angles_by_cycle_exp{i},2);
    exp_ID{i}= i*ones(curr_num_trials,1)+3;
end
close all;
all_angles_by_cycle_exp = cat(2,angles_by_cycle_exp{:});
exp_ID = cat(1,exp_ID{:});
all_swing_exp = cat(2,swing_exp{:});
all_AEP_exp = cat(2,AEP_exp{:});
all_PEP_exp = cat(2,PEP_exp{:});
N_exp = size(all_angles_by_cycle_exp,2);
x_norm = 0:0.001:1;
angle_names = {'hip'; 'knee'; 'ankle'; 'mtp'};


mean_ctr_angle = nanmean(all_angles_by_cycle_ctr,2);
std_ctr_angle = nanstd(all_angles_by_cycle_ctr,[],2);
sem_ctr_angle = std_ctr_angle/sqrt(sum(~isnan(all_angles_by_cycle_ctr(1,:,1))));
mean_swing_ctr = nanmean(all_swing_ctr);

mean_exp_angle = nanmean(all_angles_by_cycle_exp,2);
std_exp_angle = nanstd(all_angles_by_cycle_exp,[],2);
sem_exp_angle = std_exp_angle/sqrt(sum(~isnan(all_angles_by_cycle_exp(1,:,1))));
mean_swing_exp = nanmean(all_swing_exp);


for i = 1:4
    % plot mean angle throughout cycle +-sem for controls
    h(i) = figure;
    hold on
    sem_pos_ctr =  mean_ctr_angle(i,:) + sem_ctr_angle(i,:);
    sem_neg_ctr =  mean_ctr_angle(i,:) - sem_ctr_angle(i,:);
    rep_x = [x_norm, fliplr(x_norm)];
    inBetween_ctr = [sem_pos_ctr, fliplr(sem_neg_ctr)];
    fill(rep_x, inBetween_ctr,'r','FaceAlpha',.3,'EdgeAlpha',.3)
    plot(x_norm, mean_ctr_angle(i,:)','color','r');shg
    plot([mean_swing_ctr mean_swing_ctr],[0 200],'r')


   % plot mean angle throughout cycle +-sem for experimentals
    sem_pos_exp =  mean_exp_angle(i,:) + sem_exp_angle(i,:);
    sem_neg_exp =  mean_exp_angle(i,:) - sem_exp_angle(i,:);
    rep_x = [x_norm, fliplr(x_norm)];
    inBetween_exp = [sem_pos_exp, fliplr(sem_neg_exp)];
    fill(rep_x, inBetween_exp,'b','FaceAlpha',.3,'EdgeAlpha',.3)
    plot(x_norm, mean_exp_angle(i,:)','color', 'b');shg    
    % plot mean swing time (start of swing) +-std for experimentals
    plot([mean_swing_exp mean_swing_exp],[0 200],'b')
    ylabel([angle_names{i} ' angle'])
    saveas(gcf,['N:\Nofar\Behavior\DigiGait\Analysis\DTR- mean_', angle_names{i}, '_angle ctr vs exp.emf'])
end

ctr_trial_num = 1:size(all_AEP_ctr,2);
not_use_ctr = find(isnan(all_AEP_ctr(1,:)));
trial_use_ctr = setdiff(ctr_trial_num, not_use_ctr);
ctr_ID_use = ctr_ID(trial_use_ctr);
dummy_ctr = zeros(length(trial_use_ctr),1);

exp_trial_num = 1:size(all_AEP_exp,2);
not_use_exp = find(isnan(all_AEP_exp(1,:)));
trial_use_exp = setdiff(exp_trial_num, not_use_exp);
exp_ID_use = exp_ID(trial_use_exp);
dummy_exp = ones(length(trial_use_exp),1);

for i=1:4 
    all_AEP = [all_AEP_ctr(i,trial_use_ctr)'; all_AEP_exp(i,trial_use_exp)']; 
    dummy_vars = [dummy_ctr; dummy_exp];
    mouse_id = [ctr_ID_use; exp_ID_use];

    T_AEP = table(all_AEP,dummy_vars,mouse_id,'VariableNames',{'all_AEP','MouseType','MouseID'});
%     lme = fitlme(T_AEP,'all_AEP~1+MouseType+(1|MouseID)+(MouseType|MouseID)');
%     lme = fitlme(T_AEP,'all_AEP~1+MouseType+(MouseType|MouseID)')
    lme = fitlme(T_AEP,'all_AEP~1+MouseType+(1|MouseID)');
    p_AEP(i) = lme.Coefficients.pValue(2)
    clear  T* lme
    
    all_PEP = [all_PEP_ctr(i,trial_use_ctr)'; all_PEP_exp(i,trial_use_exp)']; 
    dummy_vars = [dummy_ctr; dummy_exp];
    mouse_id = [ctr_ID_use; exp_ID_use];

    T_PEP = table(all_PEP,dummy_vars,mouse_id,'VariableNames',{'all_PEP','MouseType','MouseID'});
%     lme = fitlme(T_PEP,'all_PEP~1+MouseType+(1|MouseID)+(MouseType|MouseID)');
%     lme = fitlme(T,'angle~1+MouseType+(MouseType|MouseID)')
    lme = fitlme(T_PEP,'all_PEP~1+MouseType+(1|MouseID)');
    p_PEP(i) = lme.Coefficients.pValue(2);
    clear  T* lme
end



mean_ctr_AEP = nanmean(all_AEP_ctr,2);
std_ctr_AEP = nanstd(all_AEP_ctr,[],2);
sem_ctr_AEP = nanstd(all_AEP_ctr,[],2)/sqrt(size(all_AEP_ctr,2));
mean_ctr_PEP = nanmean(all_PEP_ctr,2);
std_ctr_PEP = nanstd(all_PEP_ctr,[],2);
sem_ctr_PEP = nanstd(all_PEP_ctr,[],2)/sqrt(size(all_PEP_ctr,2));

mean_exp_AEP = nanmean(all_AEP_exp,2);
std_exp_AEP = nanstd(all_AEP_exp,[],2);
sem_exp_AEP = nanstd(all_AEP_exp,[],2)/sqrt(size(all_AEP_exp,2));
mean_exp_PEP = nanmean(all_PEP_exp,2);
std_exp_PEP = nanstd(all_PEP_exp,[],2);
sem_exp_PEP = nanstd(all_PEP_exp,[],2)/sqrt(size(all_PEP_exp,2));



bar1_loc = 0.5;
bar2_loc = 1.5;
a_ctr = bar1_loc-0.75/2;
b_ctr = bar1_loc+0.75/2;

a_exp = bar2_loc-0.75/2;
b_exp = bar2_loc+0.75/2;
for i=1:4
   % compare Anterior extreme position between ctr and exp
    figure;
    hold on
    bar(bar1_loc,mean_ctr_AEP(i),'k','BarWidth',0.75);
    bar(bar2_loc,mean_exp_AEP(i),'r','BarWidth',0.75)
    r_ctr = (b_ctr-a_ctr).*rand(N_ctr,1) + a_ctr;
    r_exp = (b_exp-a_exp).*rand(N_exp,1) + a_exp;
    plot(r_ctr, all_AEP_ctr(i,:),'o','MarkerSize',8)
    plot(r_exp, all_AEP_exp(i,:),'o','MarkerSize',8)
    errorbar(bar1_loc,mean_ctr_AEP(i),sem_ctr_AEP(i),'k');
    errorbar(bar2_loc,mean_exp_AEP(i),sem_exp_AEP(i),'k');
   
    ylabel([angle_names{i} '-AEP'])
    min_val = floor(min([all_AEP_ctr(i,:) all_AEP_exp(i,:)])/10)*10;
    max_val = ceil(max([all_AEP_ctr(i,:) all_AEP_exp(i,:)])/10)*10;
    ylim([min_val max_val])
%     xlim([-0.25 2.25])
end

for i=1:4
    % compare posterior extreme position between ctr and exp    
    figure
    hold on
    bar(bar1_loc,mean_ctr_PEP(i),'k','BarWidth',0.75);
    bar(bar2_loc,mean_exp_PEP(i),'r','BarWidth',0.75);
    errorbar(bar1_loc,mean_ctr_PEP(i),sem_ctr_PEP(i),'k');
    errorbar(bar2_loc,mean_exp_PEP(i),sem_exp_PEP(i),'k');
    r_ctr = (b_ctr-a_ctr).*rand(N_ctr,1) + a_ctr;
    r_exp = (b_exp-a_exp).*rand(N_exp,1) + a_exp;
    plot(r_ctr, all_PEP_ctr(i,:),'o','MarkerSize',8)
    plot(r_exp, all_PEP_exp(i,:),'o','MarkerSize',8)
    ylabel([angle_names{i} '-PEP']) 
    min_val = floor(min([all_PEP_ctr(i,:) all_PEP_exp(i,:)])/10)*10;
    max_val = ceil(max([all_PEP_ctr(i,:) all_PEP_exp(i,:)])/10)*10;
    ylim([min_val max_val])
end


