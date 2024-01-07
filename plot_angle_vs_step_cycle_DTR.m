function [angles_by_cycle_norm, swing_ind_norm,AEP_val,PEP_val, cycle_dur, stance_dur, swing_dur, step_num] =  plot_angle_vs_step_cycle_DTR(exp_name, animal_name)
% read DLC_Video_Info and look for video/ videos related to animal 
[~,txt,raw] =xlsread('N:\Undergrads\Mel\DLC_Video_Info.xlsx');
titles = txt(1,1:size(txt,2));
video_column = find(strcmp(titles, 'Mouse/Video Name'));
all_videos = txt(:,video_column);
video_row = find(contains(all_videos, animal_name));

 
 % load coordinate_S_updated for each video found and calculate angles 
accum_length = 0;
trans_inds(1) = 0;
 for i = 1:length(video_row)
    % load coordinates_S_updated for each video 
    project_name = raw{video_row(i),strncmp(titles, 'DLC_Project Name',16)};
    mouse_name = raw{video_row(i),strncmp(titles, 'Mouse/Video Name',16)};
    YY_joints = ['20' raw{video_row(i),strncmp(titles, 'project generation date for joints',34)}(1:2)];
    MM_joints = raw{video_row(i),strncmp(titles, 'project generation date for joints',34)}(4:5);
    DD_joints = raw{video_row(i),strncmp(titles, 'project generation date for joints',34)}(7:8);
    video_path = ['N:\Undergrads\Mel\' exp_name '\' project_name '-Mel-' YY_joints '-' MM_joints '-' DD_joints '\videos\', mouse_name];
    load([video_path '\coordinates_S_updated.mat']);
   
    % extract pixel x,y coordinates and calculate hip, knee, ankle and mtp
    % angles (angle mat is a matix n x 4 where n is the number of
    % coordinates/frames and 4 is the number of angles.
    
    curr_pixels_x = coordinates_S_updated.x_in_pixel;
    curr_pixels_y = coordinates_S_updated.y_in_pixel;    
    
    for j = 2:size(curr_pixels_x,2)-1
        curr_a = [curr_pixels_x(:,j-1)- curr_pixels_x(:,j) curr_pixels_y(:,j-1)- curr_pixels_y(:,j)];
        curr_b = [curr_pixels_x(:,j+1)- curr_pixels_x(:,j) curr_pixels_y(:,j+1)- curr_pixels_y(:,j)]; 
        for k = 1:size(curr_pixels_x,1)
           angles_cell{i}(k,j-1) = rad2deg(atan2(abs(det([curr_a(k,:);curr_b(k,:)])),dot(curr_a(k,:),curr_b(k,:)))); 
        end
    end
%   stance and swing frames of all videos from the same animal are combined
    stance_inds{i} = (coordinates_S_updated.stance_inds + accum_length)';
    swing_inds{i}= (coordinates_S_updated.swing_inds +accum_length)';
    Num_Step_Cycles(i) =  length(coordinates_S_updated.stance_inds)-1;
    max_step_dur(i)= max(diff(coordinates_S_updated.stance_inds))+1;
    accum_length =  size(curr_pixels_x,1)+accum_length;
    trans_inds(i+1)= length(stance_inds{i})+trans_inds(i);
end

% combine angles from all videos of animal
angles_mat = cat(1,angles_cell{:});
stance_inds_mat = cell2mat(stance_inds);
swing_inds_mat = cell2mat(swing_inds);
tot_Num_Step_Cycles = sum(Num_Step_Cycles);
tot_max_step_dur = max(max_step_dur);
trans_inds = trans_inds(~(trans_inds==0));

% sort angles by step cycles 
angles_by_cycle = nan(size(angles_mat,2),length(stance_inds_mat)-1,tot_max_step_dur);
percent_of_step = nan(size(angles_mat,2),length(stance_inds_mat)-1,tot_max_step_dur);

for l=1:size(angles_mat,2)
     for i= 1:length(stance_inds_mat)-1
        curr_step_cycle = stance_inds_mat(i):stance_inds_mat(i+1);
        curr_step_cycle_size = length(curr_step_cycle);
        if sum(i==trans_inds)>0 
            angles_by_cycle(l,i,:)=nan;
            swing_ind_norm(l,i) = nan;
        else
            curr_swing= swing_inds_mat(swing_inds_mat< curr_step_cycle(end)& swing_inds_mat> curr_step_cycle(1));
            angles_by_cycle(l,i,1:curr_step_cycle_size) = angles_mat(curr_step_cycle,l);    
            percent_of_step(l,i,1:curr_step_cycle_size) = 0:1/(curr_step_cycle_size-1):1;
            swing_ind_per_cycle(l,i) = curr_swing-curr_step_cycle(1)+1;
            swing_ind_norm(l,i) = percent_of_step(l,i,swing_ind_per_cycle(l,i));
        end
    end
end
% calculate the normalize transition to swing for each normalized step
% cycle
 swing_ind_norm =  swing_ind_norm(1,:);
 swing_ind_norm(swing_ind_norm==0)=nan;
 mean_swing_norm = nanmean(swing_ind_norm);
 swing_ind_per_cycle =  swing_ind_per_cycle(1,:);
 swing_ind_per_cycle(swing_ind_per_cycle==0)=nan;
 
 
% normalize each step cycle to 0-1 
% interpolate angles corresponding to x and calculate max/ min angle for
% each step cycle
x_norm = 0:0.001:1;

for i = 1:size(angles_by_cycle,1)
    for j = 1:size(angles_by_cycle,2)
          curr_size = sum(~isnan(angles_by_cycle(i,j,:)));
          if curr_size==0
              angles_by_cycle_norm(i,j,:) =nan;
          else
            for k = 1:length(x_norm)
                [~,curr_ind] = find(percent_of_step(i,j,:)==x_norm(k));
                if curr_ind
                   angles_by_cycle_norm(i,j,k) =  angles_by_cycle(i,j,curr_ind);
                else
                   angles_by_cycle_norm(i,j,k) = interp1(squeeze(percent_of_step(i,j,1:curr_size)),squeeze(angles_by_cycle(i,j,1:curr_size)),x_norm(k));
                end
            end
          end
%           angles_change_by_cycle_norm(i,j,:) = angles_by_cycle_norm(i,j,:)- angles_by_cycle_norm(i,j,1);
    end
end

%calculate AEP, PEP

for i = 1:size(angles_by_cycle,1)
    for j = 1:size(angles_by_cycle,2)
        if  isnan(swing_ind_per_cycle(j))
            PEP_val(i,j) = nan;
            AEP_val(i,j) = nan;
        else
             PEP_val(i,j) = angles_by_cycle(i,j,swing_ind_per_cycle(j));  
             AEP_val(i,j) = angles_by_cycle_norm(i,j,1001);
        end
    end
end

% calculate the mean angle throughout the step cycle
mean_angle = squeeze(nanmean(angles_by_cycle_norm,2));
std_angle = squeeze(nanstd(angles_by_cycle_norm,[],2));

% plot the mean angle (+-std) throught the step cycle for hip, knee, ankle
% and mtp 

% for i = 1:size(angles_by_cycle_norm,1)
% figure
%     std_pos =  mean_angle(i,:) + std_angle(i,:);
%     std_neg =  mean_angle(i,:) - std_angle(i,:);
%     rep_x = [x_norm, fliplr(x_norm)];
%     inBetween = [std_pos, fliplr(std_neg)];
%     fill(rep_x, inBetween,'b','FaceAlpha',.3,'EdgeAlpha',.3)
%     ylim([0 200])
%     ylabel([coordinates_S_updated.joint_names{i+1} ' angle'])
%     hold on 
%     plot(x_norm,mean_angle(i,:)','b');shg
%     line([mean_swing_norm mean_swing_norm], [-100 200],'color','b')
% %     line([mean_swing_norm+std_swing_norm mean_swing_norm+std_swing_norm], [0 200],'color','b')
% %     line([mean_swing_norm-std_swing_norm mean_swing_norm-std_swing_norm], [0 200],'color','b')
% end

% calculate the duration of swing, stance and the step cycle 
frame_rate = 415; %frames per sec
for i = 1:length(stance_inds_mat)-1
        curr_step_cycle = stance_inds_mat(i):stance_inds_mat(i+1);
        if sum(curr_step_cycle(1)==stance_inds_mat(trans_inds))>0 
            cycle_dur(i) = nan;
            swing_dur(i) = nan;
            stance_dur(i) = nan;
        else
           cycle_dur(i) = (curr_step_cycle(end)- curr_step_cycle(1)+1)/frame_rate;
           curr_swing_ind = swing_inds_mat(swing_inds_mat>curr_step_cycle(1) & swing_inds_mat<curr_step_cycle(end));
           swing_dur(i) = (curr_step_cycle(end)-  curr_swing_ind +1)/frame_rate;
           stance_dur(i) = cycle_dur(i) - swing_dur(i);    
        end
        
        step_num = size(angles_by_cycle,2);
end

   