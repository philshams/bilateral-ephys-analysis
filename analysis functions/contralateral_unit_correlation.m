% -----------------------------------------------------------------
%               Contralateral Unit Correlations
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
time_window = [-.75,1.25];
psth_bin_size = 0.001;
i = 1;



% loop across side of stimulus presentation
for stim_side = 1:2

    spike_times_timeline = ephys_data.ephys_data{stim_side,3};
    spikeDepths = ephys_data.ephys_data{stim_side,4} - (3840-insertion_depth);  
    templateDepths = ephys_data.ephys_data{stim_side,8};
    [sorted_depths,sort_ind] = sort(templateDepths);
    
    spike_clusters = ephys_data.ephys_data{stim_side,5};
    cluster_IDs  = ephys_data.ephys_data{stim_side,6};
    cluster_IDs = cluster_IDs(sort_ind);
    
    contra_stim_side = contra_sites(stim_side);
    other_spike_times_timeline = ephys_data.ephys_data{contra_stim_side ,3};
    other_spikeDepths = ephys_data.ephys_data{contra_stim_side ,4} - (3840-insertion_depth);  
    other_templateDepths = ephys_data.ephys_data{contra_stim_side ,8};
    [other_sorted_depths,other_sort_ind] = sort(other_templateDepths);
    
    other_spike_clusters = ephys_data.ephys_data{contra_stim_side ,5};
    other_cluster_IDs  = ephys_data.ephys_data{contra_stim_side ,6};
    other_cluster_IDs = other_cluster_IDs(other_sort_ind);    
        
    
    
    all_frs = cell(length(cluster_IDs),1);
    all_frs_contra = cell(length(other_cluster_IDs),1);
    
    for stim_num = 1:length(stims_by_side{stim_side})
            % align times by stimulus onset
            stim = stims_by_side{stim_side}(stim_num);
            align_times = stimOn_times(ismember(stimIDs,stim));   
            
            % loop through units
            for unit = 1:length(cluster_IDs)
                [avg_frs, ~] = get_avg_frs(spike_times_timeline(spike_clusters==cluster_IDs(unit)), ...
                                                    align_times, epoch_to_correlate, 1);
                all_frs{unit} = [all_frs{unit}; zscore(avg_frs)]; 
            end             
            
            
            for unit = 1:length(other_cluster_IDs)
                [avg_frs_contra, ~] = get_avg_frs(other_spike_times_timeline(other_spike_clusters==other_cluster_IDs(unit)), ...
                                                    align_times, epoch_to_correlate, 1);
                all_frs_contra{unit} = [all_frs_contra{unit}; zscore(avg_frs_contra)]; 
            end                           
            
    end
     
    
    % generate correlogram
    corr_coeffs_depth = zeros(length(cluster_IDs), length(other_cluster_IDs), 2);
    for unit = 1:length(cluster_IDs)
        for other_cell = 1:length(other_cluster_IDs)

            [R, P] = corrcoef(all_frs{unit},all_frs_contra{other_cell});

            if ~isfinite(R(2))
                R(2) = 0; P(2) = 1;
            end
            
            corr_coeffs_depth(unit, other_cell,1) = R(2);
            corr_coeffs_depth(unit, other_cell,2) = P(2);

        end
    end  
            
              
    
    ccp = corr_coeffs_depth(:,:,2);
    ccp = ccp(:);
    ccr = corr_coeffs_depth(:,:,1);
    ccr = ccr(:);   
    
    mean_R_contra = mean(ccr(ccp<(.01)));
    std_R_contra = std(ccr(ccp<(.01)));
    num_sig_R_contra = sum(ccp<(.01));
    frac_sig_R_contra = sum(ccp<(.01)&ccr>.1) / length(ccp);
    i = i + 1;
    
    % plot correlogram

if show_correlation_plot
        fig = figure('Position',[500 400 700 600]);

        if stim_side == 1
            corr_coeffs_depth = permute(corr_coeffs_depth,[2 1 3]);
            x_data = 'YData'; y_data = 'XData';
        else
            x_data = 'XData'; y_data = 'YData';
        end
        h = imagesc(corr_coeffs_depth(:,:,1));  
        set(h,'AlphaData',(1-corr_coeffs_depth(:,:,2)).^40)     
        title(['Pop Activity Correlogram -- ' stim_sides{stim_side} ' stim '...
            'from ' num2str(epoch_to_correlate(1)*1000) '-' num2str(epoch_to_correlate(2)*1000) 'ms'])
        colorbar
        set(gca,'Color','k')
        caxis([-.4 .8])

        set(h, x_data, sorted_depths);
        set(h, y_data, other_sorted_depths);

        xlim([0,1010])
        ylim([-40,1010])

        xlabel(['right V1 depth'])
        ylabel(['left V1 depth'])

        text(10,-20,['avg. signif. corr = ' num2str(round(mean_R_contra,2))],'color','white')
        text(400,-20,['num. signif. corr = ' num2str(sum(ccp<(.01))) ' / ' num2str(length(ccp)) ],'color','white')        

        fig.InvertHardcopy = 'off';
        saveas(fig,[save_folder 'Pop Activity Correlogram - ' stim_sides{stim_side} ' stim '...
            'from ' num2str(epoch_to_correlate(1)*1000) '-' num2str(epoch_to_correlate(2)*1000) 'ms' '.tif'])    
end
    
end

    



