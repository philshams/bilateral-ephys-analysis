% -----------------------------------------------------------------
%                Unit Correlations
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


for site = 1:2
    spike_times_timeline = ephys_data.ephys_data{site,3};
    spikeDepths = ephys_data.ephys_data{site,4} - (3840-insertion_depth);  
    templateDepths = ephys_data.ephys_data{site,8};
    [sorted_depths,sort_ind] = sort(templateDepths);
    
    spike_clusters = ephys_data.ephys_data{site,5};
    cluster_IDs  = ephys_data.ephys_data{site,6};
    cluster_IDs = cluster_IDs(sort_ind);
    
% loop across side of stimulus presentation
for stim_side = 1:2

    all_frs = cell(length(cluster_IDs),1);
    
    for stim_num = 1:length(stims_by_side{stim_side})
            % align times by stimulus onset
            stim = stims_by_side{stim_side}(stim_num);
            align_times = stimOn_times(ismember(stimIDs,stim));   
            
            % loop through depth bins
            for unit = 1:length(cluster_IDs)
                [avg_frs, ~] = get_avg_frs(spike_times_timeline(spike_clusters==cluster_IDs(unit)), ...
                                                    align_times, epoch_to_correlate, 1);
                all_frs{unit} = [all_frs{unit}; zscore(avg_frs)];
                
            end             
    end
     
    
    % generate correlogram
    corr_coeffs_depth = zeros(length(cluster_IDs), length(cluster_IDs), 2);
    for unit = 1:length(cluster_IDs)
        for other_cell = 1:length(cluster_IDs)

            [R, P] = corrcoef(all_frs{unit},all_frs{other_cell});
            if ~isfinite(R(2))
                R(2) = 0; P(2) = 1;
            end
            corr_coeffs_depth(unit, other_cell,1) = R(2);
            corr_coeffs_depth(unit, other_cell,2) = P(2);
            
        end
    end  
            
    ccp = corr_coeffs_depth(:,:,2); ccp = ccp(:);
    ccr = corr_coeffs_depth(:,:,1); ccr = ccr(:);   
    mean_R_ipsi = mean(ccr(ccp<(.01)&ccr~=1));     
    std_R_ipsi = std(ccr(ccp<(.01)&ccr~=1));     
    num_sig_R_ipsi = sum(ccp<(.01)&ccr~=1);
    frac_sig_R_ipsi = sum(ccp<(.01)&ccr~=1&ccr>.1) / sum(ccr~=1);
    i = i + 1;
    
    % plot correlogram
       
        fig = figure('Position',[300 200 700 600]);
        h = imagesc(corr_coeffs_depth(:,:,1));
        set(h,'AlphaData',(1-corr_coeffs_depth(:,:,2)).^40)     
        title(['Pop Activity Correlogram -- ' stim_sides{stim_side} ' stim / ' site_sides{site} ' V1'...
            ' from ' num2str(epoch_to_correlate(1)*1000) '-' num2str(epoch_to_correlate(2)*1000) 'ms'])
        colorbar
        set(gca,'Color','k')
        caxis([-.4 .8])

        set(h, 'XData', sorted_depths);
        set(h, 'YData', sorted_depths);
        xlim([0,1010])
        ylim([-40,1010])

        xlabel([site_sides{site} ' V1 depth'])
        ylabel([site_sides{site} ' V1 depth'])


        text(10,-20,['avg. signif. corr = ' num2str(round(mean_R_ipsi,2))],'color','white')    
        text(400,-20,['num. signif. corr = ' num2str(sum(ccp<(.01)&ccr~=1)) ' / ' num2str(sum(ccr~=1))],'color','white')    

        fig.InvertHardcopy = 'off';
        saveas(fig,[save_folder  site_sides{site} ' V1 ' 'Pop Activity Correlogram - ' stim_sides{stim_side} ' stim '...
            ' from ' num2str(epoch_to_correlate(1)*1000) '-' num2str(epoch_to_correlate(2)*1000) 'ms' '.tif'])        
    
    if stim_side == 1
        last_ccd = corr_coeffs_depth(:,:,:);
    end
    
end

if show_correlation_plot
        fig = figure('Position',[300 200 700 600]);
        stim_diff =  (corr_coeffs_depth(:,:,1) - last_ccd(:,:,1)) * ((site==2)*2-1) ;
        h = imagesc( stim_diff );
        set(h,'AlphaData',( (1-corr_coeffs_depth(:,:,2)).*(1-last_ccd(:,:,2)) ) .^10 )
        title(['Stim-Evoked Diff in Pop Activity Correlogram -- ' site_sides{site} ' V1'...
            ' from ' num2str(epoch_to_correlate(1)*1000) '-' num2str(epoch_to_correlate(2)*1000) 'ms'])
        colorbar
        set(gca,'Color','k')
        colormap_RedWhiteBlue;
        caxis([-.1 .1])

        set(h, 'XData', sorted_depths);
        set(h, 'YData', sorted_depths);
        xlim([0,1010])
        ylim([-40,1010])

        xlabel([site_sides{site} ' V1 depth'])
        ylabel([site_sides{site} ' V1 depth'])

        sig_stim_diff = stim_diff( (corr_coeffs_depth(:,:,2) < .001) | (last_ccd(:,:,2) < .001) );
        mean_diff = mean(sig_stim_diff(:));

        text(10,-20,['avg. stim-evoked diff in signif. corr. = ' num2str(round(mean_diff,2))],'color','white')          

        fig.InvertHardcopy = 'off';
        saveas(fig,[save_folder site_sides{site} ' V1 ' 'Stim-Evoked Diff in Pop Activity Correlogram - ' ...
            ' from ' num2str(epoch_to_correlate(1)*1000) '-' num2str(epoch_to_correlate(2)*1000) 'ms' '.tif'])        
end
    
end




