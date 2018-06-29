% -----------------------------------------------------------------
%                Percent Variance by Unit
% -----------------------------------------------------------------



% set up experiment parameters
contra_sites = [2 1];
window_duration = .1;


baseline_variance = cell(4,2);

for site = sites

    f_scatter = figure('Position',[800 200 800 800]);  hold on
    percent_explained_PSTH_and_contra_all_times = [];
    percent_explained_PSTH_all_times = [];
    percent_explained_PSTH_pop_times = [];
    percent_explained_PSTH_and_contra_pop_times = [];
        
    
    
    for start_time = -2:.1:-.1
        
        time_window = [start_time start_time+window_duration];
        
    
    disp([site_sides{site} ' V1'])
    percent_explained_PSTH_all = [];
    percent_explained_PSTH_and_contra_all = [];
    
    % load data
    spike_times_timeline = ephys_data.ephys_data{site,3};
    spike_clusters = ephys_data.ephys_data{site,5};
    cluster_IDs  = ephys_data.ephys_data{site,6};
    templateDepths = ephys_data.ephys_data{site,8};
    num_neurons = length(cluster_IDs);

    % load other probe's data
    contra_site = contra_sites(site);
    other_spike_times_timeline = ephys_data.ephys_data{contra_site,3};
    other_spike_clusters = ephys_data.ephys_data{contra_site,5};
    other_cluster_IDs  = ephys_data.ephys_data{contra_site,6};
    other_num_neurons = length(other_cluster_IDs);

    % do variance calculation for population data as well
    stim = stims_by_side{site};
    align_times_all = stimOn_times(ismember(stimIDs,stim)); 
    
    [avg_frs, ~] = get_avg_frs(spike_times_timeline, align_times_all, time_window, num_neurons);
    total_variance_pop = sum( var(avg_frs) );     
       
    psth_all_stim_mean_sub_pop = zeros(length(align_times_all),1); 
    psth_all_stim_mean_sub_contra_sub_pop = zeros(length(align_times_all),1);    
    
    % loop across units 
    for unit_num = 1:num_neurons

        if unit_num > 1 && ~analyze_single_units
            break
        end
        
        curr_cluster = cluster_IDs(unit_num);
        
        % get total variance
        % align times of contralateral stimulus by stimulus onset
        stim = stims_by_side{site};
        align_times_all = stimOn_times(ismember(stimIDs,stim));
           
        [avg_frs, ~] = get_avg_frs(spike_times_timeline(spike_clusters==curr_cluster), align_times_all, time_window, 1);
        total_variance = sum( var(avg_frs) ); 
                          
        % loop across stimuli
        psth_all_stim_mean_sub = zeros(length(align_times_all),1); 
        psth_all_stim_mean_sub_contra_sub = zeros(length(align_times_all),1);
        trial_num = 1;
        
        for stim_num = 1:length(stims_by_side{site})

            % align times of contralateral stimulus by stimulus onset
            stim = stims_by_side{site}(stim_num);
            align_times = stimOn_times(ismember(stimIDs,stim));  

            % get variance explained by PSTH, for 0-500 ms
            [avg_frs, ~] = get_avg_frs(spike_times_timeline(spike_clusters==curr_cluster), align_times, time_window, 1);     
            psth_all_stim_mean_sub(trial_num: trial_num + length(align_times) - 1) = avg_frs - mean(avg_frs);
                                    
            
            if unit_num==1 % do for population as well
                
                % get variance explained by PSTH, for 0-500 ms
                [avg_frs, ~] = get_avg_frs(spike_times_timeline, align_times, time_window, num_neurons);
                psth_all_stim_mean_sub_pop(trial_num: trial_num + length(align_times) - 1) = avg_frs - mean(avg_frs);    
                
            end
            
            trial_num = trial_num + length(align_times);
        end
        
            % get variance explained by PSTH and bilateral spontaneous activity
                % do linear regression 
                    % get contralateral psth for regression
            [avg_frs_contra, ~] = get_avg_frs(other_spike_times_timeline, align_times_all, time_window - shift_back_time, other_num_neurons);
            avg_frs_contra_z = zscore(avg_frs_contra);
            
            % get linear weight at each time point and subtract it
            
%             X = avg_frs_contra; % time points by trials, contralateral population fr
            shuffle_ind = randperm(length(avg_frs_contra_z));
            X = avg_frs_contra_z(shuffle_ind);
            Y = psth_all_stim_mean_sub; % time points by trials, mean-subtracted psth

            W= inv(X' * X) * X' * Y;            
            
            psth_all_stim_mean_sub_contra_sub = Y - X * W;    
            
 
        
        if unit_num==1 %do PSTH variance for population activity
            
            % get variance explained by PSTH and bilateral spontaneous activity
                % do linear regression 
                    % get contralateral psth for regression
            shuffle_ind = randperm(length(avg_frs_contra_z));
            X = avg_frs_contra_z(shuffle_ind);            
            Y = psth_all_stim_mean_sub_pop; % time points by trials, mean-subtracted psth

            W= inv(X' * X) * X' * Y;
            psth_all_stim_mean_sub_contra_sub_pop = Y - X * W;               
            
            
            % get variance explained by PSTH alone
            psth_variance_pop = sum( var(psth_all_stim_mean_sub_pop) );
            % ...and by contralateral activity
            psth_and_contra_variance_pop = sum( var(psth_all_stim_mean_sub_contra_sub_pop) );

            percent_explained_PSTH_pop = 100* (total_variance_pop - psth_variance_pop) / total_variance_pop;
            percent_explained_PSTH_and_contra_pop = 100* (total_variance_pop - psth_and_contra_variance_pop) / total_variance_pop;

            disp(['POPULATION  -- ' num2str(percent_explained_PSTH_pop) '% explained by PSTH'...
                            ' and ' num2str(percent_explained_PSTH_and_contra_pop) '% explained by PSTH and contralateral activity'])       
        end
        
        % get variance explained by PSTH alone for single unit
        psth_variance = sum( var(psth_all_stim_mean_sub) );
        % ...and by contralateral activity
        psth_and_contra_variance = sum( var(psth_all_stim_mean_sub_contra_sub) );
        
        percent_explained_PSTH = 100* (total_variance - psth_variance) / total_variance;
        percent_explained_PSTH_and_contra = 100* (total_variance - psth_and_contra_variance) / total_variance;
        
        disp(['unit ' num2str(unit_num) ' -- ' num2str(percent_explained_PSTH) '% explained by PSTH'...
                        ' and ' num2str(percent_explained_PSTH_and_contra) '% explained by PSTH and contralateral activity'])

        percent_explained_PSTH_all = [percent_explained_PSTH_all percent_explained_PSTH];
        percent_explained_PSTH_and_contra_all = [percent_explained_PSTH_and_contra_all percent_explained_PSTH_and_contra];        
        
        
    end
    
    figure(f_scatter);
    disp(['AVG UNIT -- ' num2str(nanmean(percent_explained_PSTH_all)) '% explained by PSTH'...
                    ' and ' num2str(nanmean(percent_explained_PSTH_and_contra_all)) '% explained by PSTH and contralateral activity'])    
    
    scatter(percent_explained_PSTH_all, percent_explained_PSTH_and_contra_all, 15, [.4 .3 1], 'filled',...
                            'MarkerEdgeColor', [.7 .7 .7],'MarkerFaceAlpha',.6)
    scatter(percent_explained_PSTH_pop, percent_explained_PSTH_and_contra_pop, 200, [.3 .7 1], 'filled',...
                            'MarkerEdgeColor', [.7 .7 .7],'MarkerFaceAlpha',.4)
    set(gca,'Color','k')
    xlabel('% variance explained by PSTH')
    ylabel('% variance explained by PSTH and contralateral activity')            
    title(['variance explained -- ' site_sides{site} ' V1'...
             ' from ' num2str(time_window(1)*1000) '-' num2str(time_window(2)*1000) 'ms']);
    y_limit = ylim;

    plot([xlim], [xlim], 'linestyle', ':', 'LineWidth', .5, 'color', [.9 .9 .9]);

    pause(.1);
    
    percent_explained_PSTH_all_times = [percent_explained_PSTH_all_times; percent_explained_PSTH_all];
    percent_explained_PSTH_and_contra_all_times = [percent_explained_PSTH_and_contra_all_times; percent_explained_PSTH_and_contra_all];
    percent_explained_PSTH_pop_times = [percent_explained_PSTH_pop_times percent_explained_PSTH_pop];
    percent_explained_PSTH_and_contra_pop_times = [percent_explained_PSTH_and_contra_pop_times percent_explained_PSTH_and_contra_pop];
    
    end
    
    baseline_variance{1,site} = percent_explained_PSTH_all_times;
    baseline_variance{2,site} = percent_explained_PSTH_and_contra_all_times;
    baseline_variance{3,site} = percent_explained_PSTH_pop_times;
    baseline_variance{4,site} = percent_explained_PSTH_and_contra_pop_times;
    
end


% save('baseline_variance','baseline_variance')





