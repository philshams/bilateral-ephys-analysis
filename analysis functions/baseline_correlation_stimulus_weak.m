% -----------------------------------------------------------------
%                Get Baseline Correlation Stats
% -----------------------------------------------------------------

% for every unit, get and save the coefficient, correlation, and
% significance of that unit's baseline activity with that of the contra population

% set up experiment parameters


time_window = [0 .5];
contra_sites = [2 1];
psth_bin_size = .001;

baseline_correlations_stim_weak = cell(2,2);

for site = sites

%     f = figure('Position',[800 200 800 800]);  hold on
    disp([site_sides{site} ' V1'])
    
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
   
    baseline_correlations_stim_weak{1,site} = zeros(num_neurons,3);
    baseline_correlations_stim_weak{1,site} = zeros(1,3);    
    
        % align times of contralateral stimulus by stimulus onset
        stim = stims_by_side{site};
        align_times_all = stimOn_times(ismember(stimIDs,stim));
        
        
% do for population as well    
% get variance explained by stimulus, for 0-500 ms
[pop_psth, psth_smooth, ~, ~, ~,~,~] = ... 
psth_and_smooth(spike_times_timeline, align_times_all, time_window, psth_bin_size, smooth_window, num_neurons);
            
% get contralateral psth for regression
[contra_psth, psth_smooth, ~, ~, ~,~,~] = ... 
    psth_and_smooth(other_spike_times_timeline, align_times_all, time_window-shift_back_time, psth_bin_size, smooth_window, other_num_neurons);


X =  zscore(contra_psth(:)); %-mean(contra_psth(:));% time points by trials, contralateral population fr
Y = pop_psth(:) - mean(pop_psth(:)); % time points by trials, mean-subtracted frs

[R, P] = corrcoef(X,Y);
W= inv(X' * X) * X' * Y;         

baseline_correlations_stim_weak{2,site}(1) = W;    
baseline_correlations_stim_weak{2,site}(2) = R(2);    
baseline_correlations_stim_weak{2,site}(3) = P(2);  

disp(R(2));


    % loop across units 
    for unit_num = 1: num_neurons*analyze_single_units

        curr_cluster = cluster_IDs(unit_num);

        % get psth for all stim
        [unit_psth, psth_smooth, ~, ~, ~,~,~] = ... 
            psth_and_smooth(spike_times_timeline(spike_clusters==curr_cluster), align_times_all, time_window, psth_bin_size, smooth_window, 1);
                

            
        % get variance explained by stimulus and bilateral spontaneous activity                
        X = zscore(contra_psth(:)); %-mean(contra_psth(:)); % time points by trials, contralateral population fr
        Y = unit_psth(:) - mean(unit_psth(:)); % time points by trials, mean-subtracted frs
        
        [R, P] = corrcoef(X,Y);
        W= inv(X' * X) * X' * Y;
        
        baseline_correlations_stim_weak{1,site}(unit_num,1) = W;    
        baseline_correlations_stim_weak{1,site}(unit_num,2) = R(2);   
        baseline_correlations_stim_weak{1,site}(unit_num,3) = P(2);    
        
        disp(R(2));

    end

%     histogram(baseline_correlations{1,site}(:,2))
    
end

% mean(baseline_correlations{1,1}(:,2))
% mean(baseline_correlations{1,2}(:,2))
% std(baseline_correlations{1,1}(:,2))
% std(baseline_correlations{1,2}(:,2))

% 
% for site = 1:2
%     
%     
% figure; hold on
% 
% 
% sig_baseline_corr = baseline_correlations{1,site}(baseline_correlations{1,site}(:,3)<.001,2);
% sig_stim_corr = baseline_correlations{1,site}(baseline_correlations{1,site}(:,3)<.001,2);
% 
% insig_baseline_corr = baseline_correlations{1,site}(baseline_correlations{1,site}(:,3)>=.01,2);
% insig_stim_corr = baseline_correlations{1,site}(baseline_correlations{1,site}(:,3)>=.01,2);
% 
% 
% scatter(sig_baseline_corr,sig_stim_corr, 25, [0 0 1], 'filled', 'MarkerEdgeColor', [.7 .7 .7])
% scatter(insig_baseline_corr,insig_stim_corr, 25, [.1 .1 .5], 'filled', 'MarkerEdgeColor', [1 1 1])
% 
% 
% ylim([-.1 .3])
% xlim([-.1 .3])
% 
% line(ylim,ylim,'linestyle','--','color',[.7 .7 1]);
% line([0 0],ylim,'linestyle','--','color',[.3 .3 .6]);
% line(xlim,[0 0],'linestyle','--','color',[.3 .3 .6]);
% 
% xlabel('unit to contra population correlation -- baseline')
% ylabel('unit to contra population correlation -- stimulus')
% title([site_sides{site} ' V1 unit to contra population correlations, during baseline vs stimulus'])
% set(gca,'Color','k')
% 
% end
% 
% for site = 1:2
%     
%     
% figure; hold on
% 
% sig_baseline_corr = baseline_correlations{1,site}(baseline_correlations{1,site}(:,3)<.001,2);
% 
% insig_baseline_corr = baseline_correlations{1,site}(baseline_correlations{1,site}(:,3)>=.01,2);
% 
% sig_baseline_corr_ipsi = baseline_correlations_ipsi{1,site}(baseline_correlations{1,site}(:,3)<.001,2);
% 
% insig_baseline_corr_ipsi = baseline_correlations_ipsi{1,site}(baseline_correlations{1,site}(:,3)>=.01,2);
% 
% 
% 
% scatter(sig_baseline_corr_ipsi,sig_baseline_corr, 25, [0 0 1], 'filled', 'MarkerEdgeColor', [.7 .7 .7])
% scatter(insig_baseline_corr_ipsi,insig_baseline_corr, 25, [.1 .1 .5], 'filled', 'MarkerEdgeColor', [1 1 1])
% 
% 
% ylim([-.1 .3])
% xlim([-.1 .3])
% 
% line(ylim,ylim,'linestyle','--','color',[.7 .7 1]);
% line([0 0],ylim,'linestyle','--','color',[.3 .3 .6]);
% line(xlim,[0 0],'linestyle','--','color',[.3 .3 .6]);
% 
% xlabel('unit to ipsi population correlation -- baseline')
% ylabel('unit to contra population correlation -- baseline')
% title([site_sides{site} ' V1 unit to contra population correlations, during baseline vs stimulus'])
% set(gca,'Color','k')
% 
% end
