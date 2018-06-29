% -----------------------------------------------------------------
%                Percent Variance by Unit
% -----------------------------------------------------------------


            
% set up experiment parameters
contra_sites = [2 1];


for site = sites

    f = figure('Position',[800 200 800 800]);  hold on
    disp([site_sides{site} ' V1'])
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
    
    [avg_frs, ~] = get_avg_frs(spike_times_timeline, align_times_all, epoch_to_split_by, num_neurons);
    total_variance_pop = sum( var(avg_frs) );     
       
    frs_all_stim_mean_sub_pop = zeros(size(avg_frs)); 
    frs_all_stim_mean_sub_contra_sub_pop = zeros(size(avg_frs));    
    
    % loop across units 
    for unit_num = 1:num_neurons

        curr_cluster = cluster_IDs(unit_num);
        
        % get total variance
        % align times of contralateral stimulus by stimulus onset
        stim = stims_by_side{site};
        align_times_all = stimOn_times(ismember(stimIDs,stim));
           
        [avg_frs, ~] = get_avg_frs(spike_times_timeline(spike_clusters==curr_cluster), align_times_all, epoch_to_split_by, 1);
        total_variance = sum( var(avg_frs) ); 
                          
        % loop across stimuli
        frs_all_stim_mean_sub = zeros(length(align_times),1); 
        contra_frs_all_stim = zeros(length(align_times),1); 
        frs_all_stim_mean_sub_contra_sub = zeros(length(align_times),1);
        stims_used = zeros(length(align_times),1);
        trial_num = 1;
        
        for stim_num = 1:length(stims_by_side{site})

            % align times of contralateral stimulus by stimulus onset
            stim = stims_by_side{site}(stim_num);
            align_times = stimOn_times(ismember(stimIDs,stim));  

            % get frs for each stim
            [avg_frs, ~] = get_avg_frs(spike_times_timeline(spike_clusters==curr_cluster), align_times, epoch_to_split_by, 1);     
            frs_all_stim_mean_sub(trial_num: trial_num + length(align_times) - 1) = avg_frs - mean(avg_frs);           
            stims_used(trial_num: trial_num + length(align_times) - 1) = stim; 
            
            % get contralateral frs for regression
            [avg_frs_contra, ~] = get_avg_frs(other_spike_times_timeline, align_times, epoch_to_split_by - shift_back_time, other_num_neurons);
            contra_frs_all_stim(trial_num: trial_num + length(align_times) - 1) = zscore(avg_frs_contra);
            

            if unit_num==1 % do for population as well
                % get variance explained by stimulus, for 0-500 ms
                [avg_frs, ~] = get_avg_frs(spike_times_timeline, align_times, epoch_to_split_by, num_neurons);
                frs_all_stim_mean_sub_pop(trial_num: trial_num + length(align_times) - 1) = avg_frs - mean(avg_frs);    
            end
            
            trial_num = trial_num + length(align_times);
        end
        
            % get variance explained by stimulus and bilateral spontaneous activity
                % do linear regression 
              % split between different stimulus sides   

              
            if variance_by_baseline
                X = contra_frs_all_stim; % time points by trials, contralateral population fr
                Y = frs_all_stim_mean_sub; % time points by trials, mean-subtracted frs
                W= baseline_correlations{1,site}(unit_num,1);          
                frs_all_stim_mean_sub_contra_sub = Y - X * W;    
                
                
            else % do by all trials together
                X = contra_frs_all_stim; % time points by trials, contralateral population fr
                Y = frs_all_stim_mean_sub; % time points by trials, mean-subtracted frs
                W= inv(X' * X) * X' * Y;            
                frs_all_stim_mean_sub_contra_sub = Y - X * W;    
            end
            
 
        
        if unit_num==1 %do frs variance for population activity
            
            
            % get variance explained by stimulus and bilateral spontaneous activity
                % do linear regression 
                    % get contralateral frs for regression
            if variance_by_baseline 
                X = contra_frs_all_stim; % time points by trials, contralateral population fr
                Y = frs_all_stim_mean_sub_pop; % time points by trials, mean-subtracted frs
                W= baseline_correlations{2,site}(1);          
                frs_all_stim_mean_sub_contra_sub_pop = Y - X * W;                    
                
            else % do by all trials together
                X = contra_frs_all_stim; % time points by trials, contralateral population fr
                Y = frs_all_stim_mean_sub_pop; % time points by trials, mean-subtracted frs
                W= inv(X' * X) * X' * Y;            
                frs_all_stim_mean_sub_contra_sub_pop = Y - X * W;    
            end
            
            
            
            
            
            % get variance explained by stimulus alone
            frs_variance_pop = sum( var(frs_all_stim_mean_sub_pop) );
            % ...and by contralateral activity
            frs_and_contra_variance_pop = sum( var(frs_all_stim_mean_sub_contra_sub_pop) );

            percent_explained_frs_pop = 100* (total_variance_pop - frs_variance_pop) / total_variance_pop - mean(baseline_variance{3,site});
            if variance_by_baseline
            percent_explained_frs_and_contra_pop = 100* (total_variance_pop - frs_and_contra_variance_pop) / total_variance_pop - mean(baseline_variance{3,site});                
            else
            percent_explained_frs_and_contra_pop = 100* (total_variance_pop - frs_and_contra_variance_pop) / total_variance_pop - mean(baseline_variance{4,site});                
            end


            disp(['POPULATION  -- ' num2str(percent_explained_frs_pop) '% explained by stimulus'...
                            ' and ' num2str(percent_explained_frs_and_contra_pop) '% explained by stimulus and contralateral activity'])       
        end
        
        % get variance explained by stimulus alone for single unit
        frs_variance = sum( var(frs_all_stim_mean_sub) );
        % ...and by contralateral activity
        frs_and_contra_variance = sum( var(frs_all_stim_mean_sub_contra_sub) );
        
        
        percent_explained_frs = 100* (total_variance - frs_variance) / total_variance - nanmean(baseline_variance{1,site}(:,unit_num));
        if variance_by_baseline
        percent_explained_frs_and_contra = 100* (total_variance - frs_and_contra_variance) / total_variance - nanmean(baseline_variance{1,site}(:,unit_num));    
        else
        percent_explained_frs_and_contra = 100* (total_variance - frs_and_contra_variance) / total_variance - nanmean(baseline_variance{2,site}(:,unit_num));
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
    
    scatter(percent_explained_frs_all, percent_explained_frs_and_contra_all, 15, [.4 .3 1], 'filled',...
                            'MarkerEdgeColor', [.7 .7 .7],'MarkerFaceAlpha',.6)
    scatter(percent_explained_frs_pop, percent_explained_frs_and_contra_pop, 200, [.3 .7 1], 'filled',...
                            'MarkerEdgeColor', [.7 .7 .7],'MarkerFaceAlpha',.4)
    set(gca,'Color','k')
    xlabel('% variance explained by stimulus')
    ylabel('% variance explained by stimulus and contralateral activity')            
    title(['variance explained -- ' site_sides{site} ' V1'...
             ' from ' num2str(epoch_to_split_by(1)*1000) '-' num2str(epoch_to_split_by(2)*1000) 'ms']);

    xlim([-10 80]);
    ylim([-10 80]);
	
    y_limit = ylim;
    plot([xlim], [xlim], 'linestyle', ':', 'LineWidth', .5, 'color', [.9 .9 .9]);
    
    x_limit = xlim;
    x_start = x_limit / 12;
    text(y_limit(1) + 5, y_limit(2) - (y_limit(2)-y_limit(1))/40,['-- AVG UNIT -- '],'color','white')
    text(y_limit(1) + 5, y_limit(2) - 2*(y_limit(2)-y_limit(1))/40,[num2str(round(nanmean(percent_explained_frs_all),1)) '% explained by stimulus'],'color','white')
    text(y_limit(1) + 5, y_limit(2) - 3*(y_limit(2)-y_limit(1))/40,[num2str(round(nanmean(percent_explained_frs_and_contra_all),1)) '% explained by stimulus and contralateral activity'],'color','white')
 	
    text(y_limit(1) + 5, y_limit(2) - 5*(y_limit(2)-y_limit(1))/40,[' -- POPULATION  -- ' ],'color','white')
    text(y_limit(1) + 5, y_limit(2) - 6*(y_limit(2)-y_limit(1))/40,[num2str(round(percent_explained_frs_pop,1)) '% explained by stimulus'],'color','white')
    text(y_limit(1) + 5, y_limit(2) - 7*(y_limit(2)-y_limit(1))/40,[num2str(round(percent_explained_frs_and_contra_pop,1)) '% explained by stimulus and contralateral activity'],'color','white')



    
f = figure('Position',[680 365 893 733]); hold on

set(gca,'Color','k');

% stimulus p-value vs correlation to contra
scatter(percent_explained_frs_all,percent_explained_frs_and_contra_all-percent_explained_frs_all, 25, [0 0 1], 'filled', 'MarkerEdgeColor', [.7 .7 .7])
xlabel('% variance explained by stimulus')
ylabel('% additional variance explained by contralateral side')
title([site_sides{site} ' V1 correlation between stimulus and contralateral activity explanatory power'])
% and plot best fit line
coefficients = polyfit(percent_explained_frs_all, percent_explained_frs_and_contra_all-percent_explained_frs_all, 1);
xFit = linspace(min(percent_explained_frs_all), max(percent_explained_frs_all), 100);
yFit = polyval(coefficients , xFit);

[R,P] = corrcoef(percent_explained_frs_all, percent_explained_frs_and_contra_all-percent_explained_frs_all);

hold on;
plot(xFit, yFit, 'linestyle', '--', 'LineWidth', 2, 'color', [.5 .5 .5]);
x_limit = xlim; y_limit = ylim;
text(x_limit(2) - x_limit(2) / 5,y_limit(2) - (y_limit(2)-y_limit(1))/20,['r = ' num2str(round(R(2),2))],'color','white')
text(x_limit(2) - x_limit(2) / 5,y_limit(2) - (y_limit(2)-y_limit(1))/10,['p = ' num2str(round(P(2),2))],'color','white')    


    
pause(.1);
    

end









