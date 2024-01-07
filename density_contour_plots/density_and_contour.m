
% experiment = 'PVFlpO;Lbx1Cre;Ai65' %or 'PVTdTomato'; 
% addpath(genpath('V:\Nofar\Rostrocudal analysis\final spreadsheets for rostrocaudal analysis'))
function  density_and_contour(experiment)

    segments = {'UC' 'LC' 'MT' 'UL' 'LL' 'S'};
    
    for i = 1:length(segments)
        segment = segments{i};
        coord_data = xlsread(['V:\Nofar\Rostrocudal analysis\final spreadsheets for rostrocaudal analysis\' experiment '\' segment '\' segment ' spots.csv']);
        outline_data = xlsread(['V:\Nofar\Rostrocudal analysis\final spreadsheets for rostrocaudal analysis\' experiment '\' segment '\' segment ' outline.csv']);

        h{i} = figure;
        plot(outline_data(:,1),outline_data(:,2),'o')


        [f,xi] = ksdensity(coord_data(:,1:2));
        length_xy_matrix = sqrt(length(f));
        xx = reshape(xi(:,1),length_xy_matrix ,length_xy_matrix );
        yy = reshape(xi(:,2),length_xy_matrix ,length_xy_matrix );
        zz = reshape(f(:),length_xy_matrix ,length_xy_matrix );
        normalization_factor = median(diff(unique(xi(:,1))))*median(diff(xi(:,2)));
        zz_perc = zz*normalization_factor*100;

        hold on
        contourf(xx,yy,zz_perc,6,'LineColor','none');colorbar
        xlim([-0.2 1]); 
        ylim([0 1.2]);
        yticks([0:0.1:1.2])
        title([experiment '- ' segment]);shg
        c_limits = caxis;
        caxis_min(i) = c_limits(1);
        caxis_max(i) = c_limits(2);
    end
%             saveas(h,['V:\Nofar\Rostrocudal analysis\final spreadsheets for rostrocaudal analysis\' experiment '\' segment '\' experiment '- ' segment '.emf'])
    caxis_min_all = floor(min(caxis_min)*100)/100;
    caxis_max_all = ceil(max(caxis_max)*100)/100;
    if strcmp(experiment ,'PVFlpO;Lbx1Cre;Ai65')
        map = multigradient([1 1 1 ;1 0 0]);
    else
        map = multigradient([1 1 1 ;0 0 1]);
    end
    
    for i=1:length(h)
        figure(h{i})
        colormap(map)
        caxis([caxis_min_all caxis_max_all])
        colorbar('XTick',caxis_min_all:0.25:caxis_max_all);
        saveas(h{i},['V:\Nofar\Rostrocudal analysis\final spreadsheets for rostrocaudal analysis\' experiment '\' segments{i} '\' experiment '- ' segments{i} '.emf'])
    end
end
