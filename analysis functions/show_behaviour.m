% -----------------------------------------------------------------
%                Behaviour across Session
% -----------------------------------------------------------------



% set up experiment parameters
stims_by_side = {[2,3,1,4,5], [7,8,6,9,10]};
contrasts = {'12%','25%','50%','100%'};
stim_colors = {[.45 0 .45],[0 .7 0]}; %left right
tertile_color = {[.1 .1 .8],[.5 .5 .8],[1 1 1 .9]};
contra_sites = [2 1];
trial_num = (1:length(stimOn_times)) * (stimOn_times(end) / length(stimOn_times));

% loop over stimulus etc. and plot behaviour raster

bfig = figure('Name','Behaviour Plor','Position',[144 470 1723 243]); hold on
set(gca,'Color','k') 

for stim_side = 1:2


    
site = stim_side; % site 1 (right probe) corresponds to stim_side 1 (left stimulus)
    
% load data
spike_times_timeline = ephys_data.ephys_data{site,3};
spike_clusters = ephys_data.ephys_data{site,5};
cluster_IDs  = ephys_data.ephys_data{site,6};
templateDepths = ephys_data.ephys_data{site,8};
num_neurons = length(cluster_IDs);

% load other probe's data if separating by contralateral activity
contra_site = contra_sites(site);
other_spike_times_timeline = ephys_data.ephys_data{contra_site,3};
other_spike_clusters = ephys_data.ephys_data{contra_site,5};
other_cluster_IDs  = ephys_data.ephys_data{contra_site,6};
other_num_neurons = length(other_cluster_IDs);
    
    
    for stim_num = 1:length(stims_by_side{stim_side})
    
        stim = stims_by_side{stim_side}(stim_num);
        
%         curr_stim_times = stimOn_times(ismember(stimIDs,stim))';
        curr_stim_times = trial_num(ismember(stimIDs,stim))';
        
        if max(epoch_to_split_by) > 0
            activity_stim_times = curr_stim_times;
        else
            activity_stim_times = trial_num(ismember(stimIDs,1:10))';
        end
            
        curr_raster_times = repelem(curr_stim_times, 3);
        curr_raster_times(3:3:end) = NaN;
        
        curr_raster_y = ones(size(curr_raster_times));
        curr_raster_y(2:2:end) = curr_raster_y(2:2:end) + 1;

        h_b{stim_side} = plot(curr_raster_times / 60,curr_raster_y, 'linewidth',3, 'color', [stim_colors{stim_side} .5+stim_num/10]);
        
        % get activity tertile
        [avg_frs_ipsi, trial_tertiles_ipsi] = get_avg_frs(other_spike_times_timeline, activity_stim_times, epoch_to_split_by, other_num_neurons);  
        [avg_frs_contra, trial_tertiles_contra] = get_avg_frs(spike_times_timeline, activity_stim_times, epoch_to_split_by, num_neurons);
        
        for activity_tertile = 1:3

            curr_tertile_times = activity_stim_times(trial_tertiles_contra==activity_tertile);
            curr_raster_times = repelem(curr_tertile_times, 3);
            curr_raster_times(3:3:end) = NaN;
            
            curr_raster_y = ones(size(curr_raster_times)) - site;
            curr_raster_y(2:2:end) = curr_raster_y(2:2:end) + 1;
            
            h_a{activity_tertile} = plot(curr_raster_times / 60, curr_raster_y, ...
                'linewidth',3, 'color', tertile_color{activity_tertile});
            
        end
        
        for activity_tertile = 1:3

            curr_tertile_times = activity_stim_times(trial_tertiles_ipsi==activity_tertile);
            curr_raster_times = repelem(curr_tertile_times, 3);
            curr_raster_times(3:3:end) = NaN;
            
            curr_raster_y = ones(size(curr_raster_times)) - contra_site;
            curr_raster_y(2:2:end) = curr_raster_y(2:2:end) + 1;
            
            h_a{activity_tertile} = plot(curr_raster_times / 60, curr_raster_y, ...
                'linewidth',3, 'color', tertile_color{activity_tertile});
            
        end        
        
        
        
    end
    
end

% title([stim_sides{stim_side} ' stim trial type and activity from ' num2str(epoch_to_split_by(1)*1000) '-' num2str(epoch_to_split_by(2)*1000) 'ms']);
title(['trial type and activity from ' num2str(epoch_to_split_by(1)*1000) '-' num2str(epoch_to_split_by(2)*1000) 'ms']);
set(gca,'ytick',[-.5 .5 1.5])
set(gca,'yticklabel',{'left activity' 'right activity' 'stimulus'})
xlabel('time (minutes)')
xlim([0 22])
legend([h_b{1} h_b{2} h_a{1} h_a{2} h_a{3}], {'left','right','low activity','medium activity','high activity'},'TextColor','w','location','northeastoutside')
% legend([h_a{1} h_a{2} h_a{3}], {'low activity','medium activity','high activity'},'TextColor','w','location','northeastoutside')
    
