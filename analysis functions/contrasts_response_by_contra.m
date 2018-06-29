% -----------------------------------------------------------------
%            Get Contrast Response Curve by Contra Activity
% -----------------------------------------------------------------


%% PLOT PSTH BY STIMULUS AND OPTIONALLY BY CONTRALATERAL ACTIVITY TERTILE

% set up experiment parameters
stims_by_side = {[5,4,3,2,7,8,9,10], [5,4,3,2,7,8,9,10]};

contrast_labels = {'100% L','50% L','25% L','12% L','12% R','25% R','50% R','100% R'};
contrasts_for_model = [NaN -12 -25 -50 -100 NaN 12 25 50 100];


stim_sides = {'left','right'};
tertile_color = {[.3 .3 1],[.5 .5 .85],[.7 .7 .7]};
tertile_color_2 = {[0 0 1],[.5 .5 .85],[1 1 1]};

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
 fig = figure('Name',['site ' num2str(site) ' PSTH'],'Position',[700 50 1000 1050]); hold on;
for site = sites

    stim_side = site;
    
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
    
    subplot(2,1,site); 
    
    
    
    for stim_num = 1:length(stims_by_side{stim_side})
        
        % align times by stimulus onset
        stim = stims_by_side{stim_side}(stim_num);
        align_times = stimOn_times(ismember(stimIDs,stim));   

        [frs, ~] = get_avg_frs(spike_times_timeline, align_times, time_window, num_neurons);
        [baseline_frs, ~] = get_avg_frs(spike_times_timeline, align_times, [-1 0], num_neurons);   
        stim_responses = frs - baseline_frs;

        % get trial nums of trial in current contralateral activity tertile            
        [~, trial_tertiles_contra] = get_avg_frs(other_spike_times_timeline, align_times, time_window, other_num_neurons);

        % concatenate data
        all_stim_responses{stim_num} = stim_responses;
        all_trial_tertiles_contra{stim_num} = trial_tertiles_contra;
        
        
        
    end
    

    % loop across bins of magnitude of contralateral activity
    for contralateral_activity_tertile = [1 2 3] %1:3^(separate_contralateral_spontaneous)

        contrasts = contrasts_for_model(stims_by_side{site})';
        
        curr_stim_responses = zeros(length(stims_by_side{stim_side}),2);
        contrast_response_for_model = zeros(0,2);
        
%         %do so for all trials
%         for stim_num = 1:length(stims_by_side{stim_side})
%             curr_stim_responses(stim_num,1) = mean(all_stim_responses{stim_num}(:));
%             curr_stim_responses(stim_num,2) = std(all_stim_responses{stim_num}(:)) ...
%                                                                          / sqrt(length(all_stim_responses{stim_num}(:)));
%         
%             stim = stims_by_side{stim_side}(stim_num);
%             stim_responses = all_stim_responses{stim_num}(:);
%             contrast_response_for_model = [contrast_response_for_model; [ones(size(stim_responses))*contrasts_for_model(stim) stim_responses] ];
%         end        
%         if site == 1 % set other side stimulus to zero contrast
%             contrast_response_for_model(contrast_response_for_model(:,1) > 0,1) = 0;
%         else
%             contrast_response_for_model(contrast_response_for_model(:,1) < 0,1) = 0;
%         end
%         
%         % 	[ err, pars ] = fit_hyper_ratio(cs,resps,nn,Rmax, sigma, n, R0) imposes each paramters if not empty
%         [ err, pars ] = fit_hyper_ratio_DS(abs(contrast_response_for_model(:,1)),contrast_response_for_model(:,2));
%         sigma = pars(2); n = pars(3);
        
        
        for stim_num = 1:length(stims_by_side{stim_side})
            curr_stim_responses(stim_num,1) = mean(all_stim_responses{stim_num}(all_trial_tertiles_contra{stim_num}==contralateral_activity_tertile));
            curr_stim_responses(stim_num,2) = std(all_stim_responses{stim_num}(all_trial_tertiles_contra{stim_num}==contralateral_activity_tertile)) ...
                                                                         / sqrt(length(all_stim_responses{stim_num}(all_trial_tertiles_contra{stim_num}==contralateral_activity_tertile)));
        
            stim = stims_by_side{stim_side}(stim_num);
            stim_responses = all_stim_responses{stim_num}(all_trial_tertiles_contra{stim_num}==contralateral_activity_tertile);
            contrast_response_for_model = [contrast_response_for_model; [ones(size(stim_responses))*contrasts_for_model(stim) stim_responses] ];
        end
        % plot output
%         fill([contrasts; flip(contrasts)], [(curr_stim_responses(:,1)-curr_stim_responses(:,2)); flip(curr_stim_responses(:,1)+curr_stim_responses(:,2))], ...
%             tertile_color{contralateral_activity_tertile}, 'EdgeAlpha', 0, 'FaceAlpha', .3); 
        hold on
        h{contralateral_activity_tertile} = errorbar(contrasts+contralateral_activity_tertile, curr_stim_responses(:,1), ...
                                curr_stim_responses(:,2), 'color',tertile_color{contralateral_activity_tertile}, 'linewidth',3,'marker','.','markers',20);

        % do model

        if site == 1 % set other side stimulus to zero contrast
            contrast_response_for_model(contrast_response_for_model(:,1) > 0,1) = 0;
        else
            contrast_response_for_model(contrast_response_for_model(:,1) < 0,1) = 0;
        end
        
        % 	[ err, pars ] = fit_hyper_ratio(cs,resps,nn,Rmax, sigma, n, R0) imposes each paramters if not empty
        [ err, pars ] = fit_hyper_ratio_DS(abs(contrast_response_for_model(:,1)),contrast_response_for_model(:,2),[],[],[],[],[]);
        disp(pars);
        
        x_contrast = linspace(min(contrast_response_for_model(:,1)),max(contrast_response_for_model(:,1)));        
        model_output = hyper_ratio(pars,abs(x_contrast));
        
        plot(x_contrast, model_output , 'color',tertile_color_2{contralateral_activity_tertile}, 'linewidth',2, 'linestyle',':');
        plot([0 -100*sign(x_contrast(2))],[min(model_output) min(model_output)], 'color',tertile_color_2{contralateral_activity_tertile}, 'linewidth',2, 'linestyle',':');
        
                                                            
        title(['Contrast Response Function of ' site_sides{site} ' V1 (' num2str(1000*time_window(1)) ' - ' num2str(1000*time_window(2)) 'ms)']);

        set(gca,'Color','k')
        
        ylabel('avg firing rate - baseline');
        xlabel('stimulus')
        xticks([-100 -50 -25 -12 12 25 50 100])
        xticklabels(contrast_labels)
        
    end
    locations = {'northeast','northwest'};
    legend([h{3},h{2},h{1}],{'high','medium','low'},'location',locations{site},'TextColor','w')
    xlim([-110 110]);

end