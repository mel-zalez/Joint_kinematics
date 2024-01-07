% this function calculates coordinates (pixels and mm) of the knee. 
% function input:
% exp_name-  'PV_Project_DTR' or 'PV_Project_ReachR'
% num_of_frame- when adjusting knee location,  choose a small # of frames (50-100). 
% choose 1 when done when done adjusting (function will go over all frames)
% when wanting he function to go through all frames (after figuring out
% knee location)
%save_outputs: 0 when still sdjusting knee location (video with new skeleton is not
%saved) 1- when done adjusting knee location, will generate and save a
%video with knew skeleton that includes the knee. 

 
function generate_new_skeleton(exp_name,video_name,length_f,length_t, num_of_frames, save_outputs)
    
[status,sheets] = xlsfinfo('V:\Undergrads\Mel\DLC_Video_Info.xlsx');
sheet_num = find(strcmp(exp_name,sheets));
[num,txt,raw] =xlsread('V:\Undergrads\Mel\DLC_Video_Info.xlsx',sheet_num);
titles = txt(1,1:size(txt,2));

video_column = find(strcmp(titles, 'Mouse/Video Name'));
all_videos = txt(:,video_column);
video_row = find(strcmp(all_videos, video_name));



%     find all info for joints excel sheet name
project_name = raw{video_row,strncmp(titles, 'DLC_Project Name',16)};
mouse_name = raw{video_row,strncmp(titles, 'Mouse/Video Name',16)};
YY_joints = ['20' raw{video_row,strncmp(titles, 'project generation date for joints',34)}(1:2)];
MM_joints = raw{video_row,strncmp(titles, 'project generation date for joints',34)}(4:5);
DD_joints = raw{video_row,strncmp(titles, 'project generation date for joints',34)}(7:8);
month_name_joints = datestr(datetime(1,str2num(MM_joints),1),'mmm');
iterations_joints = num2str(raw{video_row,strncmp(titles, 'joints Training Iterations',26)});
joints_file_name = ['V:\Undergrads\Mel\' exp_name '\' project_name '-Mel-' YY_joints '-' MM_joints '-' DD_joints '\videos\' mouse_name 'DLC_resnet50_' project_name month_name_joints  DD_joints 'shuffle1_' iterations_joints];


%     find all info for joints excel sheet name
YY_cal = ['20' raw{video_row,find(strncmp(titles, 'project generation date for calibration points',34))}(1:2)];
MM_cal = raw{video_row,find(strncmp(titles, 'project generation date for calibration points',34))}(4:5);
DD_cal = raw{video_row,find(strncmp(titles, 'project generation date for calibration points',34))}(7:8);
month_name_cal = datestr(datetime(1,str2num(MM_cal),1),'mmm');
iterations_cal = num2str(raw{video_row,find(strcmp(titles, 'Calibration Training Iterations '))'});
calibration_file_name = ['V:\Undergrads\Mel\' exp_name '\pixelstomm-Mel-' YY_cal '-' MM_cal '-' DD_cal '\videos\' mouse_name 'DLC_resnet50_pixelstomm' month_name_cal DD_cal 'shuffle1_' iterations_cal];

    
% read joint and calibration excel files 
[num_joints,txt_joints,raw_joints] =xlsread([joints_file_name, '.csv']);
[num_cal,txt_cal,raw_cal] =xlsread([calibration_file_name, '.csv']);
joint_names = txt_joints(2,1:size(txt_joints,2));


% find coversion to mm for each frame

x_cal = [num_cal(:,2) num_cal(:,5) num_cal(:,8) num_cal(:,11) num_cal(:,14) num_cal(:,17)];
y_cal = [num_cal(:,3) num_cal(:,6) num_cal(:,9) num_cal(:,12) num_cal(:,15) num_cal(:,18)];

diff_x = [x_cal(:,2)-x_cal(:,1) x_cal(:,4)-x_cal(:,3) x_cal(:,6)-x_cal(:,5)];
diff_y = [y_cal(:,2)-y_cal(:,1) y_cal(:,4)-y_cal(:,3) y_cal(:,6)-y_cal(:,5)];
dist = sqrt(power(diff_x,2) + power(diff_y,2));
mean_dist = mean(dist,2);
mm_in_px = mean_dist./5

length_f_px = length_f*mm_in_px;
length_t_px = length_t*mm_in_px;


% extract x (first column) and y coordiantes (second column), in pixels, for each frame for every joint.
iliac_crest_pos = find(strncmp(joint_names, 'Illiaccrest',10));
iliac_crest_px = num_joints(:,iliac_crest_pos(1:2));

hip_pos = find(strncmp(joint_names, 'Hip',3));
hip_px = num_joints(:,hip_pos(1:2));

ankle_pos = find(strncmp(joint_names, 'Ankle',5));
ankle_px = num_joints(:,ankle_pos(1:2));

mtp_pos = find(strncmp(joint_names, 'MTP',3));
mtp_px = num_joints(:,mtp_pos(1:2));

toe_pos = find(strncmp(joint_names, 'Toe',5));
toe_px = num_joints(:,toe_pos(1:2));

% calculate corrected knee location based on known distance of femur (length_f) and tibia (length_t) 
knee_px = nan(length(iliac_crest_px),2);
if num_of_frames ==1
    num_of_frames = length(iliac_crest_px);
end


for i=1:num_of_frames
    x_h = hip_px(i,1);
    y_h = hip_px(i,2);
    x_a = ankle_px(i,1);
    y_a = ankle_px(i,2);
    syms x_k y_k
    E = [(x_k - x_a)^2 + (y_k - y_a)^2 == length_t_px(i)^2,(x_k - x_h)^2 + (y_k - y_h)^2 == length_f_px(i)^2 ];
    S = solve(E,x_k,y_k);
    if vpa(S.y_k(2))>y_h && vpa(S.x_k(2))>x_a
        knee_px(i,1) = vpa(S.x_k(2));
        knee_px(i,2) = vpa(S.y_k(2));
    else 
        knee_px(i,1) = vpa(S.x_k(1));
        knee_px(i,2) = vpa(S.y_k(1));
    end
    disp(i)
    clear S E x_h y_h x_a y_a x_k y_k
end


% Access video file
v = VideoReader([joints_file_name '_labeled.mp4']);


% Preallocate structur'e to store video frames
% Note that this is just an initial guess for preallocation based on
% duration and framerate, the video may have fewer or more frames

nFrames = length(iliac_crest_px);
s(nFrames) = struct('cdata',[],'colormap',[]);

figure
colormap gray

k=1;
while hasFrame(v)
    im = readFrame(v);
    mymovie(:,:,k) = im(:,:,1);
    k=k+1;
end

for i = 1:num_of_frames
    colormap gray
    image(mymovie(:,:,i));shg
    hold on
    plot([iliac_crest_px(i,1),hip_px(i,1),knee_px(i,1),ankle_px(i,1),mtp_px(i,1),toe_px(i,1)],[iliac_crest_px(i,2),hip_px(i,2),knee_px(i,2),ankle_px(i,2),mtp_px(i,2),toe_px(i,2)],'Marker','o','MarkerSize',2,'color','r');     
    title(i)
    curframe = getframe(gcf);
    mymovie2(:,:,i) = curframe.cdata(:,:,1);
    title(i)
    hold off
    if ~save_outputs
     pause
    end
end

if save_outputs
    clear vOut
    video_path = ['V:\Undergrads\Mel\' exp_name '\' project_name '-Mel-' YY_joints '-' MM_joints '-' DD_joints '\videos\', mouse_name];
    mkdir(video_path)
    vOut = VideoWriter([video_path,'\', mouse_name] ,'MPEG-4');

    frame = getframe;
    frame.colormap = [];


    open(vOut);
    for l=1:size(mymovie2,3)
        frame.cdata = repmat(mymovie2(:,:,l),[1 1 3]);
        writeVideo(vOut,frame);
    end
    close(vOut)

    coordinates_S. joint_names = {'iliac crest', 'hip', 'knee', 'ankle', 'mtp', 'toe'};
    coordinates_S.x_in_mm = [iliac_crest_px(:,1)./mm_in_px, hip_px(:,1)./mm_in_px, knee_px(:,1)./mm_in_px, ankle_px(:,1)./mm_in_px, mtp_px(:,1)./mm_in_px, toe_px(:,1)./mm_in_px];
    coordinates_S.y_in_mm = [iliac_crest_px(:,2)./mm_in_px, hip_px(:,2)./mm_in_px, knee_px(:,2)./mm_in_px, ankle_px(:,2)./mm_in_px, mtp_px(:,2)./mm_in_px, toe_px(:,2)./mm_in_px];
    coordinates_S.x_in_pixel = [iliac_crest_px(:,1), hip_px(:,1), knee_px(:,1), ankle_px(:,1), mtp_px(:,1), toe_px(:,1)];
    coordinates_S.y_in_pixel = [iliac_crest_px(:,2), hip_px(:,2), knee_px(:,2), ankle_px(:,2), mtp_px(:,2), toe_px(:,2)];

    save([video_path, '\coordinates_S'],'coordinates_S')
end










