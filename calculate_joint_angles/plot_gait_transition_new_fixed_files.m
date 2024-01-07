function  attract_all = plot_gait_transition_new_fixed_files()
 
speed_vec = 20:20:80;
base_folder= 'N:\Behavior data\Data\DigiGait\PV project\PVFlpO;Lbx1Cre;DTR';
[~,~,raw] =xlsread([base_folder '\Copy of Paw area_all_animals_Mel2.0.xlsx']);

speeds = nan(length(raw),1);
condition = strings(length(raw),1);

for i=2:length(raw)
    speeds(i) = raw{i,4};
    condition(i) = raw{i,6};
    exclude(i,1) = raw{i,10};
    fixed{i,1} = raw{i,12};
end

condition_types = ['C' , 'E'];

for i=1:length(condition_types)
   for j = 1:length(speed_vec)
    curr_speed = speed_vec(j);
    curr_ind = find(condition==condition_types(i) & speeds==curr_speed & exclude==0 & strcmp(fixed, 'fixed')); 
        for m = 1:length(curr_ind)
                curr_dob= raw{curr_ind(m),3};
                slash_ind = find(curr_dob=='/');
                animal_name = raw{curr_ind(m),1};
                file_name = [base_folder, '\PVFlpO;Lbx1Cre;DTR;Ai65 ',curr_dob(1:slash_ind(1)-1),'-',curr_dob(slash_ind(1)+1:slash_ind(2)-1),'-20\Paw contact area new\PawArea_',animal_name,'_',num2str(curr_speed),'cms_0degUP'];
                curr_num = xlsread(file_name);
                [~,~,~,~,~,~,~,~,~,occur_mat(m,:), persist_mat(m,:), attract_mat_temp] = get_circular_limb_coupling_values(curr_num);
                attract_mat(m,:,:) = attract_mat_temp;
                clear curr_num attract_mat_temp
        end
          occur_all{i,j} = nanmean(occur_mat);
          persist_all{i,j} = nanmean(persist_mat);
          attract_all{i,j} = squeeze(nanmean(attract_mat,1));
          clear occur_mat persist_mat attract_mat
   end
end

gaits = {'Walk'; 'Trot';'Gallop';'Bound'};
gait_attrac_indx = [2 3 4 ; 1 3 4 ; 1 2 4 ; 1 2 3]'; 



rgb = [1 1 1 ; 0 0 0];
map = multigradient(rgb,'length',100);


for i=1:length(condition_types)
   for j = 1:length(speed_vec)
       
        %define arrow width 2 for 100% occurence
        occur_norm = 0.02;
        curr_occur = occur_all{i,j}*occur_norm;
        curr_persist = persist_all{i,j};
        figure 
        hold on
        pos_W = [0 0 curr_occur(1) curr_occur(1)];
        rectangle('Position',pos_W,'Curvature',[1 1],'FaceColor', map(ceil(curr_persist(1))+1,:));        
        pos_Trot = [4 0 curr_occur(2) curr_occur(2)];
        rectangle('Position',pos_Trot,'Curvature',[1 1],'FaceColor', map(ceil(curr_persist(2))+1,:));
        pos_Gallop = [4 4 curr_occur(3) curr_occur(3)];
        rectangle('Position',pos_Gallop,'Curvature',[1 1],'FaceColor', map(ceil(curr_persist(3))+1,:));
        pos_Bound = [0 4 curr_occur(4) curr_occur(4)];
        rectangle('Position',pos_Bound,'Curvature',[1 1],'FaceColor', map(ceil(curr_persist(4))+1,:));
        xlim([-1 6])
        ylim([-1 5])
        axis equal;shg
        
        %define arrow width 30 for 100% attraction
        attrac_norm = 0.3;
%         draw arrows to LW
       
        if attract_all{i,j}(1,1)>0
           arrow([3.5 0.3],[1.5 0.3],'Width',attract_all{i,j}(1,1)*attrac_norm);%from Trot to Walk
        end
        if attract_all{i,j}(2,1)>0
           arrow([3.5 3.5],[1.5,1.5],'Width',attract_all{i,j}(2,1)*attrac_norm);%from Gallop to Walk
        end
        if attract_all{i,j}(3,1)>0
           arrow([0.5 3.5],[0.5 2],'Width',attract_all{i,j}(3,1)*attrac_norm);%from Bound to Walk
        end

        
%         draw arrows to Trot
        if attract_all{i,j}(1,2)>0
           arrow([1.5 0.4],[3.5 0.4],'Width',attract_all{i,j}(1,2)*attrac_norm);%from Walk to Trot
        end
        if attract_all{i,j}(2,2)>0
           arrow([4.2 3.5],[4.2 1.5],'Width',attract_all{i,j}(2,2)*attrac_norm);%from Gallop to Trot
        end      
        if attract_all{i,j}(3,2)>0
           arrow([1.5 3.5],[3.5 1.5],'Width',attract_all{i,j}(3,2)*attrac_norm);%from Bound to Trot
        end 
          
        %         draw arrows to Gallop
        if attract_all{i,j}(1,3)>0
           arrow([1.4 1.6],[3.4 3.6],'Width',attract_all{i,j}(1,3)*attrac_norm);%from Walk to Gallop
        end
        
        if attract_all{i,j}(2,3)>0
           arrow([4.1 1.5],[4.1 3.5],'Width',attract_all{i,j}(2,3)*attrac_norm);%from Trot to Gallop
        end
        
        if attract_all{i,j}(3,3)>0
           arrow([1.5 4],[3.5 4],'Width',attract_all{i,j}(3,3)*attrac_norm);%from Bound to Gallop
        end
        
        %         draw arrows to Bound
        if attract_all{i,j}(1,4)>0
           arrow([0.6 2],[0.6 3.5],'Width',attract_all{i,j}(1,4)*attrac_norm);%from Walk to Bound
        end
        if attract_all{i,j}(2,4)>0
           arrow([3.6 1.6],[1.6 3.6],'Width',attract_all{i,j}(2,4)*attrac_norm);%from Trot to Bound
        end

        if attract_all{i,j}(3,4)>0
           arrow([3.5 3.9],[1.5 3.9],'Width',attract_all{i,j}(3,4)*attrac_norm);%from Gallop to Bound
        end  
        
        colormap(map)
        colorbar;shg
        if i==1
            title(['ctr_ at_ ' num2str(speed_vec(j)) ' _cm/s'])
%             saveas(gcf,['V:\Nofar\Behavior\DigiGait\Analysis\ctr_gait_trans_speed' num2str(speed_vec(j)) '.emf'])

        else
            title(['exp_ at_ ' num2str(speed_vec(j)) ' _cm/s'])
%             saveas(gcf,['V:\Nofar\Behavior\DigiGait\Analysis\exp_gait_trans_speed' num2str(speed_vec(j)) '.emf'])
        end       
   end
end

figure
pos_W = [0 0 occur_norm*100 occur_norm*100];
rectangle('Position',pos_W,'Curvature',[1 1])
xlim([-1 6])
ylim([-1 5])
axis equal;shg
% saveas(gcf,'V:\Nofar\Behavior\DigiGait\Analysis\100_per_occur.emf')

      
figure
xlim([-1 6])
ylim([-1 5])
axis equal;shg
arrow([2 1.5],[2 3.5], 'Width',attrac_norm*100)
% saveas(gcf,'V:\Nofar\Behavior\DigiGait\Analysis\100_per_attract.emf')
