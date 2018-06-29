% -----------------------------------------------------------------
%            Get Contrast Response Curve by Contra Activity
% -----------------------------------------------------------------


%% PLOT PSTH BY STIMULUS AND OPTIONALLY BY CONTRALATERAL ACTIVITY TERTILE

% set up experiment parameters
stims_by_side = {[5,4,3,2,7,8,9,10], [5,4,3,2,7,8,9,10]};
contrast_labels = {'100% L','50% L','25% L','12% L','12% R','25% R','50% R','100% R'};
stim_sides = {'left','right'};
tertile_color = {[.2 .2 .6],[.4 .4 .8],[.7 .7 1]};
contra_sites = [2 1];
avg_frs_contra_all = [];
avg_frs_all = [];
frs_contra_all = [];
frs_all = [];
h = {};
time_window = epoch_to_split_by;
psth_bin_size = 0.001;
stimulus_colors = [];
corr_fig = {}; hm = {};

% loop over probes
 fig = figure('Name',['site ' num2str(site) ' PSTH'],'Position',[700 50 800 1050]); hold on;
for site = sites

    all_stim_responses = cell(length(stims_by_side{stim_side}),1);
    all_trial_tertiles_contra = cell(length(stims_by_side{stim_side}),1);

    % load data
    spike_times_timeline = ephys_data.ephys_data{site,3};
    cluster_IDs  = ephys_data.ephys_data{site,6};
    num_neurons = length(cluster_IDs);

    % load other probe's data if separating by contralateral activity
    contra_site = contra_sites(site);
    other_spike_times_timeline = ephys_data.ephys_data{contra_site,3};
    other_cluster_IDs  = ephys_data.ephys_data{contra_site,6};
    other_num_neurons = length(other_cluster_IDs);
    contra_text = [' by contra activity during ' num2str(epoch_to_split_by(1)*1000) '-' num2str(epoch_to_split_by(2)*1000) 'ms'];

    % use contralateral stimulus 
    stim_side = site;
    subplot(2,1,site); 
    
    
    
    for stim_num = 1:length(stims_by_side{stim_side})
        
        % align times by stimulus onset
        stim = stims_by_side{stim_side}(stim_num);
        align_times = stimOn_times(ismember(stimIDs,stim));   

        [frs, ~] = get_avg_frs(spike_times_timeline, align_times, time_window, num_neurons);
        [baseline_frs, ~] = get_avg_frs(spike_times_timeline, align_times, [-1 0], num_neurons);   
        stim_responses = frs - baseline_frs;

        % concatenate data
        all_stim_responses{stim_num} = stim_responses;
    end

    

    % loop across bins of magnitude of contralateral activity
    contrasts = (1:length(stims_by_side{stim_side}))';

    curr_stim_responses = zeros(length(stims_by_side{stim_side}),2);
    for stim_num = 1:length(stims_by_side{stim_side})
        curr_stim_responses(stim_num,1) = mean(all_stim_responses{stim_num});
        curr_stim_responses(stim_num,2) = std(all_stim_responses{stim_num});
    end
    % plot output
    fill([contrasts; flip(contrasts)], [(curr_stim_responses(:,1)-curr_stim_responses(:,2)); flip(curr_stim_responses(:,1)+curr_stim_responses(:,2))], ...
        tertile_color{3}, 'EdgeAlpha', 0, 'FaceAlpha', .3); 
    hold on
    plot(contrasts, curr_stim_responses(:,1), 'color',tertile_color{3}, 'linewidth',2,...
                                                            'marker','.','markers',20);

    title(['Contrast Response Function of ' site_sides{site} ' V1 (' num2str(1000*time_window(1)) ' - ' num2str(1000*time_window(2)) 'ms)']);

    set(gca,'Color','k')

    ylabel('avg firing rate');
    xlabel('stimulus')
    xticklabels(contrast_labels)
        


end