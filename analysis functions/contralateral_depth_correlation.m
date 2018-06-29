% -----------------------------------------------------------------
%                Contralateral Depth Correlations
% -----------------------------------------------------------------

%% PLOT CORRELOGRAM OF ACTIVITY BY DEPTH

% set up experiment parameters

tertile_color = {[.2 .2 .6],[.4 .4 .8],[.7 .7 1]};
contra_sites = [2 1];
avg_frs_contra_all = [];
avg_frs_all = [];
frs_contra_all = [];
frs_all = [];
h = {};

psth_bin_size = 0.001;
depths = 0:100:1000;




    
    
% loop across side of stimulus presentation
for stim_side = 1:2

    spike_times_timeline = ephys_data.ephys_data{stim_side,3}; % gives site contralateral to stimulus
    spikeDepths = ephys_data.ephys_data{stim_side,4} - (3840-insertion_depth);  
    cluster_IDs  = ephys_data.ephys_data{stim_side,6};
    templateDepths = ephys_data.ephys_data{stim_side,8};
    num_neurons = length(cluster_IDs);

    contra_stim_side = contra_sites(stim_side);
    other_spike_times_timeline = ephys_data.ephys_data{contra_stim_side,3};
    other_spikeDepths = ephys_data.ephys_data{contra_stim_side,4} - (3840-insertion_depth); 
    other_cluster_IDs  = ephys_data.ephys_data{contra_stim_side,6};
    other_templateDepths = ephys_data.ephys_data{contra_stim_side,8};
    other_num_neurons = length(other_cluster_IDs);        
    
    
    for stim_num = 1:length(stims_by_side{stim_side})
        
        
    all_frs = cell(length(depths)-1,1);
    all_frs_contra = cell(length(depths)-1,1);
    
    for stim_num = 1:length(stims_by_side{stim_side})
            % align times by stimulus onset
            stim = stims_by_side{stim_side}(stim_num);
            align_times = stimOn_times(ismember(stimIDs,stim));   
            
            % loop through depth bins
            for other_depth_bin = 1:(length(depths)-1)
                num_depth_neurons_other_depth = sum(other_templateDepths>=depths(other_depth_bin)&other_templateDepths<=depths(other_depth_bin+1));

                [avg_frs_contra, ~] = get_avg_frs(other_spike_times_timeline(other_spikeDepths>=depths(other_depth_bin)&...
                                                                    other_spikeDepths<=depths(other_depth_bin+1)), ...
                                                                    align_times, epoch_to_correlate, num_depth_neurons_other_depth);  

                all_frs_contra{other_depth_bin} = [all_frs_contra{other_depth_bin}; zscore(avg_frs_contra)];
            end
            
            for depth_bin = 1:(length(depths)-1)
                num_depth_neurons = sum(templateDepths>=depths(depth_bin)&templateDepths<=depths(depth_bin+1));
                [avg_frs, ~] = get_avg_frs(spike_times_timeline(spikeDepths>=depths(depth_bin)&spikeDepths<=depths(depth_bin+1)), ...
                                                    align_times, epoch_to_correlate, num_depth_neurons );
                all_frs{depth_bin} = [all_frs{depth_bin}; zscore(avg_frs)]; 
            end             
    end
    
    end
    
    % generate correlogram
    corr_coeffs_depth = zeros(length(depths)-1, length(depths)-1, 2);
    for depth_bin = 1:(length(depths)-1)
        for other_depth_bin = 1:(length(depths)-1)
            num_depth_neurons = sum(templateDepths>=depths(depth_bin)&templateDepths<=depths(depth_bin+1));
            num_depth_neurons_contra = sum(other_templateDepths>=depths(other_depth_bin)&other_templateDepths<=depths(other_depth_bin+1));
           
            if num_depth_neurons && num_depth_neurons_contra
                [R, P] = corrcoef(all_frs{depth_bin},all_frs_contra{other_depth_bin});
                corr_coeffs_depth(depth_bin, other_depth_bin,1) = R(2);
                corr_coeffs_depth(depth_bin, other_depth_bin,2) = P(2);
            else
                corr_coeffs_depth(depth_bin, other_depth_bin,1) = 0;
                corr_coeffs_depth(depth_bin, other_depth_bin,2) = 1;
            end
        end
    end  
            
                
                         
    figure('Position',[200 100 700 600]);
    h = imagesc(corr_coeffs_depth(:,:,1));
	set(h,'AlphaData',(1-corr_coeffs_depth(:,:,2)).^10)     
    title(['Pop Activity Correlogram -- ' stim_sides{stim_side} ' stim '...
        'from ' num2str(epoch_to_correlate(1)*1000) '-' num2str(epoch_to_correlate(2)*1000) 'ms'])
    caxis([-.4 .7]);
    colorbar
    set(gca,'Color','k')
    
    set(h, 'XData', depths(2:end) - depths(2)/2);
    set(h, 'YData', depths(2:end) - depths(2)/2);
    xlim([depths(1),depths(end)])
    ylim([depths(1),depths(end)])
    
    ylabel([site_sides{stim_side} ' V1 depth'])
    xlabel([stim_sides{stim_side} ' V1 depth'])
    

end

    
