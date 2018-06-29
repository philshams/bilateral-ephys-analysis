% -----------------------------------------------------------------
%                Percent Variance by Unit
% -----------------------------------------------------------------


   

% set up experiment parameters
start_times = -.5:.025:.9;



for site = sites
    
    stim_side = site;
    pop_var_frs = [];
    pop_var_frs_contra = [];

    avg_unit_var_frs = [];
    std_unit_var_frs = [];
    avg_unit_var_frs_contra = [];
    std_unit_var_frs_contra = [];
    
    med_unit_var_frs = [];
    iqr_unit_var_frs = [];
    med_unit_var_frs_contra = [];
    iqr_unit_var_frs_contra = [];    
    
for time_start = start_times

    disp([site_sides{site} ' V1']);
    time_window = [time_start time_start+window_duration]
    percent_explained_frs_all = [];
    percent_explained_frs_and_contra_all = [];
        
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
    frs_all_stim_mean_sub_pop = zeros(length(avg_frs),1); 
    frs_all_stim_mean_sub_contra_sub_pop = zeros(length(avg_frs),1);    
    
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
        frs_all_stim_mean_sub = zeros(size(avg_frs)); 
        frs_all_stim_mean_sub_contra_sub = zeros(size(avg_frs));
        contra_frs_all_stim = zeros(size(avg_frs));
        trial_num = 1;
        
        for stim_num = 1:length(stims_by_side{stim_side})
            
            % align times of contralateral stimulus by stimulus onset
            stim = stims_by_side{site}(stim_num);
            align_times = stimOn_times(ismember(stimIDs,stim));  

            % get frs for each stim
            [avg_frs, ~] = get_avg_frs(spike_times_timeline(spike_clusters==curr_cluster), align_times, time_window, 1);     
            frs_all_stim_mean_sub(trial_num: trial_num + length(align_times) - 1) = avg_frs - mean(avg_frs);           
            stims_used(trial_num: trial_num + length(align_times) - 1) = stim; 
            
            % get contralateral frs for regression
            [avg_frs_contra, ~] = get_avg_frs(other_spike_times_timeline, align_times, time_window - shift_back_time, other_num_neurons);
            contra_frs_all_stim(trial_num: trial_num + length(align_times) - 1) = zscore(avg_frs_contra); % - mean(avg_frs_contra);
            

            if unit_num==1 % do for population as well
                % get variance explained by stimulus, for 0-500 ms
                [avg_frs, ~] = get_avg_frs(spike_times_timeline, align_times, time_window, num_neurons);
                frs_all_stim_mean_sub_pop(trial_num: trial_num + length(align_times) - 1) = avg_frs - mean(avg_frs);    
            end
            trial_num = trial_num + length(align_times);
        end
        
            % get variance explained by stimulus and bilateral spontaneous activity
                % do linear regression 
              % split between different stimulus sides   
            X = contra_frs_all_stim; % time points by trials, contralateral population fr
            Y = frs_all_stim_mean_sub; % time points by trials, mean-subtracted frs
            
            if variance_by_baseline
                W = baseline_correlations{1,site}(unit_num,1);          
            else % do by all trials together
                W = inv(X' * X) * X' * Y;         
                
            end
                    
            frs_all_stim_mean_sub_contra_sub = Y - X * W; 
            
        
        if unit_num==1 %do for population activity as well
            
            % get variance explained by stimulus and bilateral spontaneous activity
                % do linear regression 
                    % get contralateral frs for regression
                X = contra_frs_all_stim; % time points by trials, contralateral population fr
                Y = frs_all_stim_mean_sub_pop; % time points by trials, mean-subtracted frs                    
            if variance_by_baseline
                W= baseline_correlations{2,site}(1);          
            else % do by all trials together
                W = inv(X' * X) * X' * Y; 
            end
             frs_all_stim_mean_sub_contra_sub_pop = Y - X * W;                     
            
            % get variance explained by stimulus alone
            stim_variance_pop = sum( var(frs_all_stim_mean_sub_pop) );
            % ...and by contralateral activity
            stim_and_contra_variance_pop = sum( var(frs_all_stim_mean_sub_contra_sub_pop) );

            percent_explained_frs_pop = 100* (total_variance_pop - stim_variance_pop) / total_variance_pop - mean(baseline_variance{3,site});
            if variance_by_baseline
            percent_explained_frs_and_contra_pop = 100* (total_variance_pop - stim_and_contra_variance_pop) / total_variance_pop - mean(baseline_variance{3,site});    
            else
            percent_explained_frs_and_contra_pop = 100* (total_variance_pop - stim_and_contra_variance_pop) / total_variance_pop - mean(baseline_variance{4,site});
            end
            disp(['POPULATION  -- ' num2str(percent_explained_frs_pop) '% explained by stimulus'...
                            ' and ' num2str(percent_explained_frs_and_contra_pop) '% explained by stimulus and contralateral activity'])       
        end
        
        % get variance explained by stimulus alone
        stim_variance = sum( var(frs_all_stim_mean_sub) );
        % ...and by contralateral activity
        stim_and_contra_variance = sum( var(frs_all_stim_mean_sub_contra_sub) );
        
        percent_explained_frs = 100* (total_variance - stim_variance) / total_variance - nanmean(baseline_variance{1,site}(:,unit_num));
        if variance_by_baseline
        percent_explained_frs_and_contra = 100* (total_variance - stim_and_contra_variance) / total_variance - nanmean(baseline_variance{1,site}(:,unit_num));
        else
        percent_explained_frs_and_contra = 100* (total_variance - stim_and_contra_variance) / total_variance - nanmean(baseline_variance{2,site}(:,unit_num));
        end
        
        disp(['unit ' num2str(unit_num) ' -- ' num2str(percent_explained_frs) '% explained by stimulus'...
                        ' and ' num2str(percent_explained_frs_and_contra) '% explained by stimulus and contralateral activity'])

        if isfinite(percent_explained_frs) && isfinite(percent_explained_frs_and_contra)
            percent_explained_frs_all = [percent_explained_frs_all percent_explained_frs];
            percent_explained_frs_and_contra_all = [percent_explained_frs_and_contra_all percent_explained_frs_and_contra];        
        end
    end

    disp(['AVG UNIT -- ' num2str(nanmean(percent_explained_frs_all)) '% explained by stimulus'...
                    ' and ' num2str(nanmean(percent_explained_frs_and_contra_all)) '% explained by stimulus and contralateral activity'])    

                
    pop_var_frs = [pop_var_frs percent_explained_frs_pop];
    pop_var_frs_contra = [pop_var_frs_contra percent_explained_frs_and_contra_pop];

    avg_unit_var_frs = [avg_unit_var_frs nanmean(percent_explained_frs_all)];
    std_unit_var_frs = [std_unit_var_frs nanstd(percent_explained_frs_all)];
    avg_unit_var_frs_contra = [avg_unit_var_frs_contra nanmean(percent_explained_frs_and_contra_all)];
    std_unit_var_frs_contra = [std_unit_var_frs_contra nanstd(percent_explained_frs_and_contra_all)];                
                
    med_unit_var_frs = [med_unit_var_frs nanmedian(percent_explained_frs_all)];
    iqr_unit_var_frs = [iqr_unit_var_frs iqr(percent_explained_frs_all)];
    med_unit_var_frs_contra = [med_unit_var_frs_contra nanmedian(percent_explained_frs_and_contra_all)];
    iqr_unit_var_frs_contra = [iqr_unit_var_frs_contra iqr(percent_explained_frs_and_contra_all)];                   
    
    
end

    time_axis = (start_times + window_duration)*1000;
    stim_colors = {[.75 0 .75],[0 .7 0]}; %left right
    
    % plot pop var
    f = figure('Position',[669 339 1068 574]); hold on
    set(gca,'Color','k') 

    p1 = plot(time_axis, pop_var_frs,'color',[.4 .3 1 .5],'linewidth',3);
    p2 = plot(time_axis, pop_var_frs_contra,'color',[.3 .7 1 .5],'linewidth',3, 'linestyle', '--');

    ylim([-10 80])
    line([0,0],ylim,'linestyle','--','color',[.7 .7 1]);
    line([500,500],ylim,'linestyle','--','color',[.5 .5 .8]);
            
    title(['population activity variance explained over time, ' site_sides{site} ' V1 ' title_suffix])
    legend([p1 p2],{'stimulus alone','stimulus and contralateral activity'},'TextColor','w')
    ylabel('% variance explained')
    xlabel('time (ms)')
    
    if analyze_single_units 
%     % plot avg unit var
    f = figure('Position',[669 339 1068 574]); hold on
    set(gca,'Color','k') 

    p1 = plot(time_axis, avg_unit_var_frs,'color',[.4 .3 1 .8],'linewidth',3);
    fill([time_axis flip(time_axis)], [(avg_unit_var_frs-std_unit_var_frs) flip(avg_unit_var_frs+std_unit_var_frs)], ...
        [.4 .3 1], 'EdgeAlpha', 0, 'FaceAlpha', .3); 
        
    p2 = plot(time_axis, avg_unit_var_frs_contra,'color',[.3 .7 1 .8],'linewidth',3);
    fill([time_axis flip(time_axis)], [(avg_unit_var_frs_contra-std_unit_var_frs_contra) flip(avg_unit_var_frs_contra+std_unit_var_frs_contra)], ...
        [.3 .7 1 ], 'EdgeAlpha', 0, 'FaceAlpha', .3); 
    
    line([0,0],ylim,'linestyle','--','color',[.7 .7 1]);
    line([500,500],ylim,'linestyle','--','color',[.5 .5 .8]);
            
    title(['avg single unit activity variance explained over time, ' site_sides{site} ' V1 ' title_suffix])
    legend([p1 p2],{'stimulus alone','stimulus and contralateral activity'},'TextColor','w')
    ylabel('% variance explained')
    xlabel('time (ms)')
    

    % plot pop var
    f = figure('Position',[669 339 1068 574]); hold on
    set(gca,'Color','k') 

    p1 = plot(time_axis, med_unit_var_frs,'color',[.4 .3 1 .8],'linewidth',3);
    fill([time_axis flip(time_axis)], [(med_unit_var_frs-iqr_unit_var_frs/2) flip(med_unit_var_frs+iqr_unit_var_frs/2)], ...
        [.4 .3 1], 'EdgeAlpha', 0, 'FaceAlpha', .3); 
        
    p2 = plot(time_axis, med_unit_var_frs_contra,'color',[.3 .7 1 .8],'linewidth',3);
    fill([time_axis flip(time_axis)], [(med_unit_var_frs_contra-iqr_unit_var_frs_contra/2) flip(med_unit_var_frs_contra+iqr_unit_var_frs_contra/2)], ...
        [.3 .7 1 ], 'EdgeAlpha', 0, 'FaceAlpha', .3); 
    
    line([0,0],ylim,'linestyle','--','color',[.7 .7 1]);
    line([500,500],ylim,'linestyle','--','color',[.5 .5 .8]);
            
    title(['avg single unit activity variance explained over time, ' site_sides{site} ' V1 ' title_suffix])
    legend([p1 p2],{'stimulus alone','stimulus and contralateral activity'},'TextColor','w')
    ylabel('% variance explained')
    xlabel('time (ms)')
    
    end
    
    pause(.1)

end










