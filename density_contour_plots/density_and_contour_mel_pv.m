% experiment = 'PVFlpO;Lbx1Cre;Ai65' %or 'PVTdTomato';
% addpath(genpath('V:\Nofar\Rostrocudal analysis\final spreadsheets for rostrocaudal analysis'))
function  density_and_contour_mel_pv(experiment)

%     segments = {'UC' 'LC' 'MT' 'UL' 'LL' 'S'};
    segments = {'LL'};
    
    for i = 1:length(segments)
        segment = segments{i};
        %coord_data = xlsread(['N:\Nofar\Rostrocudal analysis\final spreadsheets for rostrocaudal analysis\' experiment '.csv']);
        coord_data = xlsread(['V:\Nofar\Rostrocudal analysis\final spreadsheets for rostrocaudal analysis\' experiment '\' segment '\' segment ' spots.csv']);
        outline_data = xlsread(['V:\Nofar\Rostrocudal analysis\final spreadsheets for rostrocaudal analysis\' experiment '\' segment '\' segment ' outline.csv']);
        %outline_data = xlsread(['N:\Nofar\Rostrocudal analysis\final spreadsheets for rostrocaudal analysis\Overlay2' '.csv']);

        h{i} = figure;
        plot(outline_data(:,1),outline_data(:,2),'k') %,'LineWidth',2)

        [f,xi] = ksdensity(coord_data(:,1:2));
        length_xy_matrix = sqrt(length(f));
        xx = reshape(xi(:,1),length_xy_matrix ,length_xy_matrix );
        yy = reshape(xi(:,2),length_xy_matrix ,length_xy_matrix );
        zz = reshape(f(:),length_xy_matrix ,length_xy_matrix );
        normalization_factor = median(diff(unique(xi(:,1))))*median(diff(xi(:,2)));
        zz_perc = zz*normalization_factor*100;

        hold on %to include sc outline from line 14
        contour(xx,yy,zz_perc);colorbar %shaded outline of sc
%         contourf(xx,yy,zz_perc);colorbar %NO sc outline
%         xlim([-0.2 1]); 
%         ylim([0 1.2]);
        xlim([0 1]); 
        ylim([0 1]);
        yticks([0:0.2:1])%change max to 1.2
        title([experiment '-' segment]);shg
        colorbar
        c_limits = caxis;
        caxis_min(i) = c_limits(1);
        caxis_max(i) = c_limits(2);
        
   end
%             saveas(h,['V:\Nofar\Rostrocudal analysis\final spreadsheets for rostrocaudal analysis\' experiment '\' segment '\' experiment '- ' segment '.emf'])
    caxis_min_all = floor(min(caxis_min)*100)/100;
    caxis_max_all = ceil(max(caxis_max)*100)/100;
    if strcmp(experiment ,'PVFlpO;Lbx1Cre;Ai65')
        map = multigradient([1 1 1 ;1 0 0]); %red [1 0 0] for PVFlp0
    else
%         map = multigradient([1 1 1 ;.3 .3 .3]); %gray [.3 .3 .3] for TdTomatoe
%          map = multigradient([1 1 1 ; .164 .72 .46]); %green
        map = multigradient([1 1 1 ; 1 .7 .256]); %orange
    end
    
    for i=1:length(h)
        figure(h{i})
        colormap(map)
        caxis([caxis_min_all caxis_max_all])
        colorbar('XTick',caxis_min_all:0.25:caxis_max_all);
%         saveas(h{i},['V:\Nofar\Rostrocudal analysis\final spreadsheets for rostrocaudal analysis\' experiment '\' segments{i} '\' experiment '- ' segments{i} '.emf'])
    end
end