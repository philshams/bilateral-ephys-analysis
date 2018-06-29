% -----------------------------------------------------------------
%                Depth Correlations
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


for site = 1:2
    spike_times_timeline = ephys_data.ephys_data{site,3};
    spikeDepths = ephys_data.ephys_data{site,4} - (3840-insertion_depth);  
    templateDepths = ephys_data.ephys_data{site,8};
   
    
% loop across side of stimulus presentation
for stim_side = 1:2

    all_frs = cell(length(depths)-1,1);
    
    for stim_num = 1:length(stims_by_side{stim_side})
            % align times by stimulus onset
            stim = stims_by_side{stim_side}(stim_num);
            align_times = stimOn_times(ismember(stimIDs,stim));   
            
            % loop through depth bins
            for depth_bin = 1:(length(depths)-1)
                num_depth_neurons = sum(templateDepths>=depths(depth_bin)&templateDepths<=depths(depth_bin+1));
                [avg_frs, ~] = get_avg_frs(spike_times_timeline(spikeDepths>=depths(depth_bin)&spikeDepths<=depths(depth_bin+1)), ...
                                                    align_times, epoch_to_correlate, num_depth_neurons );
                all_frs{depth_bin} = [all_frs{depth_bin}; zscore(avg_frs)]; 
            end             
    end
     
    
    % generate correlogram
    corr_coeffs_depth = zeros(length(depths)-1, length(depths)-1, 2);
    for depth_bin = 1:(length(depths)-1)
        for other_depth_bin = 1:(length(depths)-1)
            num_depth_neurons = sum(templateDepths>=depths(depth_bin)&templateDepths<=depths(depth_bin+1));
            num_depth_neurons_other_depth = sum(templateDepths>=depths(other_depth_bin)&templateDepths<=depths(other_depth_bin+1));
            if num_depth_neurons && num_depth_neurons_other_depth
                [R, P] = corrcoef(all_frs{depth_bin},all_frs{other_depth_bin});
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
    title(['Pop Activity Correlogram -- ' stim_sides{stim_side} ' stim / ' site_sides{site} ' V1'...
        ' from ' num2str(epoch_to_correlate(1)*1000) '-' num2str(epoch_to_correlate(2)*1000) 'ms'])
    colorbar
    set(gca,'Color','k')
    caxis([-.2 1])
    
    set(h, 'XData', depths(2:end) - depths(2)/2);
    set(h, 'YData', depths(2:end) - depths(2)/2);
    xlim([depths(1),depths(end)])
    ylim([depths(1),depths(end)])
    
    xlabel([site_sides{site} ' V1 depth'])
    ylabel([site_sides{site} ' V1 depth'])

end

    
end





