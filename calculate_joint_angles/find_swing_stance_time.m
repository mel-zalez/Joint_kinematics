% example:
% exp_name = 'PV_Project_DTR';
% video_name = 'EN0_2.17.21_Wh_M_20cms_TR2_trimmed_10 steps';

function [coordinates_S_updated] = find_swing_stance_time(exp_name,video_name)
 
% get location of coordinate_S for the video in video_name[num,txt,raw] =xlsread('V:\Undergrads\Mel\DLC_Video_Info.xlsx');
[~,txt,raw] =xlsread('V:\Undergrads\Mel\DLC_Video_Info.xlsx');
titles = txt(1,1:size(txt,2));
video_column = find(strncmp(titles, 'Mouse/Video Name',16));
all_videos = txt(:,video_column);
video_row = find(strcmp(all_videos, video_name));
project_name = raw{video_row,strncmp(titles, 'DLC_Project Name',16)};
mouse_name = raw{video_row,strncmp(titles, 'Mouse/Video Name',16)};
YY_joints = ['20' raw{video_row,strncmp(titles, 'project generation date for joints',34)}(1:2)];
MM_joints = raw{video_row,strncmp(titles, 'project generation date for joints',34)}(4:5);
DD_joints = raw{video_row,strncmp(titles, 'project generation date for joints',34)}(7:8);

video_path = ['V:\Undergrads\Mel\' exp_name '\' project_name '-Mel-' YY_joints '-' MM_joints '-' DD_joints '\videos\', mouse_name];

% load coordinate_S for 
load([video_path '\coordinates_S.mat'])

% extract toe x and y coordinates (in pixels) from coordintes_S
toe_ind = find(strcmp(coordinates_S.joint_names  ,'toe'));
toe_x = coordinates_S.x_in_pixel(:,toe_ind);
toe_y = coordinates_S.y_in_pixel(:,toe_ind);

% plot found indexes on top of a graph showing the change in toe_y
f = figure('WindowState','maximized');
ax = axes; 
plot(toe_y,'k')

% here we identify steps to connect

orig_xlim = xlim; 
answer = inputdlg('How many steps to connect?');
if str2num(answer{1})>0
    for i=1:str2num(answer{1})
        xlim(orig_xlim)
        curr_xlim = inputdlg({'Enter minimal x limit', 'Enter maximal x limit'}, ['X limits of point #' i]) ;
        xlim([str2num(curr_xlim{1}) str2num(curr_xlim{2})])
        waitforbuttonpress;
        point_x(i) = ceil(ax.CurrentPoint(1,1));
        point_y(i) = ceil(ax.CurrentPoint(1,2));
        waitforbuttonpress;
        point_x1(i) = ceil(ax.CurrentPoint(1,1));
        point_y1(i) = ceil(ax.CurrentPoint(1,2));
        clear curr_xlim
    end
    toe_y(point_x+1:point_x1-1) = interp1([point_x point_x1],[point_y point_y1],point_x+1:point_x1-1);
    clear point*
end
% update plot with changes
hold off
plot(toe_y,'k')




%%find the times mouse transitioned to stance
[~,stance_inds] = findpeaks(toe_x); %here we find the approximate indx where mouse transitioned to stance 


% plot found indexes on top of a graph showing the change in toe_y
f = figure('WindowState','maximized');
ax = axes; 
plot(toe_y,'k')
hold on
plot(stance_inds,toe_y(stance_inds),'Marker','o','LineStyle','none','color','b');shg
title('stance times')

% here we can identify stance points we want to clear
orig_xlim = xlim; 
answer = inputdlg('How many contact points to clear?','add contact points to clear');
if str2num(answer{1})>0
    for i=1:str2num(answer{1})
        xlim(orig_xlim)
        curr_xlim = inputdlg({'Enter minimal x limit', 'Enter maximal x limit'}, ['X limits of point #' i]) ;
        xlim([str2num(curr_xlim{1}) str2num(curr_xlim{2})])
        waitforbuttonpress;
        point_x(i) = ax.CurrentPoint(1,1);
        [~,ind_to_clear] = min(abs(stance_inds-point_x(i)));
        stance_inds(ind_to_clear)=[];
        % update plot with changes
        hold off
        plot(toe_y,'k')
        hold on
        plot(stance_inds,toe_y(stance_inds),'Marker','o','LineStyle','none','color','b');shg
    end
    clear point*
end


% here we can identify additional frames on which the toe touched the
% ground
orig_xlim = xlim; 
answer = inputdlg('How many missing contact points?','add missing contact points');
if str2num(answer{1})>0
    for i=1:str2num(answer{1})
        xlim(orig_xlim)
        curr_xlim = inputdlg({'Enter minimal x limit', 'Enter maximal x limit'}, ['X limits of point #' i]) ;
        xlim([str2num(curr_xlim{1}) str2num(curr_xlim{2})])
        waitforbuttonpress;
        point_x(i) = ceil(ax.CurrentPoint(1,1));
        point_y(i) = ceil(ax.CurrentPoint(1,2));
        clear curr_xlim
    end
    stance_inds = sort([stance_inds; point_x']);
    clear point*
end
% update plot with changes
hold off
plot(toe_y,'k')
hold on
plot(stance_inds,toe_y(stance_inds),'Marker','o','LineStyle','none','color','b');shg


%for each indx we try to find the more accurate indx of transition
for i = 1:length(stance_inds)
    if (stance_inds(i)-5)<0 
        curr_inds = 1:stance_inds(i)+5;
    elseif (stance_inds(i)+5)>length(toe_y)
            curr_inds = stance_inds(i)-5:length(toe_y);
    else  
        curr_inds = stance_inds(i)-5:stance_inds(i)+5;
        flipped_curr_inds =  fliplr(curr_inds);
        curr_diff = diff((-toe_y( flipped_curr_inds )));
%         curr_max_diff = max(curr_diff);
%         stance_inds(i) = curr_inds(find(curr_diff==curr_max_diff)+1);
        
          m = find(curr_diff/max(abs(curr_diff))>.25);
          if isempty(m)
              continue
          elseif length(m)>1
            corrected_ind = m(find(diff(m)==1,1,'first'))+1;
            stance_inds(i) = flipped_curr_inds(corrected_ind);
          else
            corrected_ind = m;
            stance_inds(i) = flipped_curr_inds(corrected_ind);
          end
    end 
end

% update plot with changes
hold off
plot(toe_y,'k')
hold on
plot(stance_inds,toe_y(stance_inds),'Marker','o','LineStyle','none','color','b');shg


%%find the times mouse transitioned to swing
[~,swing_inds] = findpeaks(-toe_x); %here we find the approximate indx where mouse transitioned to stance 


% plot found indexes on top of a graph showing the change in toe_y
hold on
plot(swing_inds,toe_y(swing_inds),'Marker','o','LineStyle','none','color','r');shg


% here we can identify swing points we want to clear
orig_xlim = xlim; 
answer = inputdlg('How many swing points to clear?','add swing points to clear');
if str2num(answer{1})>0
    for i=1:str2num(answer{1})
        xlim(orig_xlim)
        curr_xlim = inputdlg({'Enter minimal x limit', 'Enter maximal x limit'}, ['X limits of point #' i]) ;
        xlim([str2num(curr_xlim{1}) str2num(curr_xlim{2})])
        waitforbuttonpress;
        point_x(i) = ax.CurrentPoint(1,1);
        [~,ind_to_clear] = min(abs(swing_inds-point_x(i)));
        swing_inds(ind_to_clear)=[];
        hold off
        plot(toe_y,'k')
        hold on
        plot(stance_inds,toe_y(stance_inds),'Marker','o','LineStyle','none','color','b');
        plot(swing_inds,toe_y(swing_inds),'Marker','o','LineStyle','none','color','r');shg
    end
    clear point*
end



% here we can identify additional frames on which the toe touched the
% ground
orig_xlim = xlim; 
answer = inputdlg('How many missing swing points?','add missing swing points');
if str2num(answer{1})>0
    for i=1:str2num(answer{1})
        xlim(orig_xlim)
        curr_xlim = inputdlg({'Enter minimal x limit', 'Enter maximal x limit'}, ['X limits of point #' i]) ;
        xlim([str2num(curr_xlim{1}) str2num(curr_xlim{2})])
        waitforbuttonpress;
        point_x(i) = ceil(ax.CurrentPoint(1,1));
        point_y(i) = ceil(ax.CurrentPoint(1,2));
        clear curr_xlim
    end
    swing_inds = sort([swing_inds; point_x']);
end
close(gcf)


% %for each indx we try to find the more accurate indx of transition
% n= min(length(swing_inds), length(stance_inds));
% for i = 1:n
%     if  swing_inds(i)-stance_inds(i)>50 || swing_inds(i)-stance_inds(i)<-50
%         x=8;
%     else
%         x=5;
%     end
%     if swing_inds(i)-x<0
%         curr_inds = 1:swing_inds(i)+x;
%     elseif (swing_inds(i)+x)>length(toe_y)
%         curr_inds = swing_inds(i)-x:length(toe_y);
%     else  
%         curr_inds = swing_inds(i)-x:swing_inds(i)+x;
%     end
%     curr_diff = diff((toe_y(curr_inds)));
%     if sum(curr_diff>1)==0
%         continue
%     else
%          m = find(curr_diff/max(abs(curr_diff))>.25);
%         if isempty(m)
%             continue
%         elseif length(m)>1
%             corrected_ind = m(find(diff(m)==1,1,'first'))+1;
%             swing_inds(i) = curr_inds(corrected_ind);
%         else 
%              corrected_ind = m;
%             swing_inds(i) = curr_inds(corrected_ind);
%         end
%     end
% end 
% 
hold off
plot(toe_y,'k')
hold on
plot(stance_inds,toe_y(stance_inds),'Marker','o','LineStyle','none','color','b');
plot(swing_inds,toe_y(swing_inds),'Marker','o','LineStyle','none','color','r');shg

coordinates_S_updated = coordinates_S;
coordinates_S_updated.swing_inds = swing_inds;
coordinates_S_updated.stance_inds = stance_inds;


save([video_path, '\coordinates_S_updated'],'coordinates_S_updated')
writetable(struct2table(coordinates_S_updated), 'Structure_Example.csv')

