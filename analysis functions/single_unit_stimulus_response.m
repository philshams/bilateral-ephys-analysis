% -----------------------------------------------------------------
%                Stimulus Responses
% -----------------------------------------------------------------

%% PLOT PSTH BY STIMULUS AND OPTIONALLY BY CONTRALATERAL ACTIVITY TERTILE

% set up experiment parameters

tertile_color = {[.2 .2 .8],[.5 .5 .8 .9],[.9 .9 1 .7]};
contra_sites = [2 1];
avg_frs_contra_all = [];
avg_frs_all = [];
frs_contra_all = [];
frs_all = [];
h = {};
time_window = [-1,1.5];
psth_bin_size = 0.001;
stimulus_colors = [];
stim_color = [0 0 .3; .3 .3 .7; .4 .4 .9; .8 .8 .8; 1 1 1];
corr_fig = {}; hm = {};

% loop over probes
for site = sites

avg_frs_contra_side = [];
avg_frs_side = [];
frs_contra_side = [];
frs_side = [];
frs_contra_side_mean = [];
frs_side_mean = [];


% load data
spike_times_timeline = ephys_data.ephys_data{site,3};
spike_clusters = ephys_data.ephys_data{site,5};
cluster_IDs  = ephys_data.ephys_data{site,6};
templateDepths = ephys_data.ephys_data{site,8};
num_neurons = length(cluster_IDs);

% load other probe's data if separating by contralateral activity
if separate_contralateral_spontaneous
    contra_site = contra_sites(site);
    other_spike_times_timeline = ephys_data.ephys_data{contra_site,3};
    other_spike_clusters = ephys_data.ephys_data{contra_site,5};
    other_cluster_IDs  = ephys_data.ephys_data{contra_site,6};
    other_num_neurons = length(other_cluster_IDs);
    contra_text = [' by contra activity during ' num2str(epoch_to_split_by(1)*1000) '-' num2str(epoch_to_split_by(2)*1000) 'ms'];
    cb = 0;
else
    contra_text = '';
    cb = 2;
end


% loop across units 
for unit_num = 1:num_neurons
        max_fr = [];
    curr_cluster = cluster_IDs(unit_num);
    % loop across side of stimulus presentation
    for stim_side = 1:2
         
        curr_figure_name = [site_sides{site} ' site - ' stim_sides{stim_side} ' stim - ' num2str(round(templateDepths(unit_num))) 'um - unit ' num2str(unit_num) ];
        psth_fig = figure('Name',curr_figure_name,'Position',[700 50 1200 1050]); hold on;

        for stim_num = 1:length(stims_by_side{stim_side})

                % align times by stimulus onset
                trial_offset = 0;
                stim = stims_by_side{stim_side}(stim_num);
                align_times_all = stimOn_times(ismember(stimIDs,stim));               

                % get trial nums of trial in current contralateral activity tertile
                if separate_contralateral_spontaneous               
                    [avg_frs_contra, trial_tertiles] = get_avg_frs(other_spike_times_timeline, align_times_all, epoch_to_split_by, other_num_neurons);  
                end            

            % loop across bins of magnitude of contralateral activity
            for contralateral_activity_tertile = 1:3^(separate_contralateral_spontaneous)

                subplot(length(stims_by_side{stim_side}),2,stim_num*2-1);
                if separate_contralateral_spontaneous
                    align_times = align_times_all(trial_tertiles==contralateral_activity_tertile);
                else
                    align_times = align_times_all;
                end

                % get PSTH
                [psth_smooth_all, psth_smooth, psth_std_trials, bins, rasterX,rasterY,spikeCounts] = ... 
                    psth_and_smooth(spike_times_timeline(spike_clusters==curr_cluster), align_times, time_window, psth_bin_size, smooth_window, 1);

                % plot output
                hold on
                h{contralateral_activity_tertile} = plot(bins * 1000, psth_smooth, 'color',tertile_color{contralateral_activity_tertile+cb}, 'linewidth',2);

                title([site_sides{site} ' V1 to ' stim_sides{stim_side} ' ' contrasts{stim_num} ' stimulus' contra_text]);

                set(gca,'Color','k')       
                ylim([min(psth_smooth)-1,max(psth_smooth)+1])
                ylabel('avg firing rate');
                xlabel('Time from stim onset (ms)')
                
                if separate_contralateral_spontaneous
                    fill([epoch_to_split_by*1000 flip(epoch_to_split_by*1000)], [-100 -100 1000 1000], ...
                    [.4 .4 .4], 'EdgeAlpha', 0, 'FaceAlpha', .07); 
                end
                
                subplot(length(stims_by_side{stim_side}),2,stim_num*2); hold on
                plot(rasterX*1000,trial_offset + rasterY,'color',tertile_color{contralateral_activity_tertile+cb},'linewidth',2)
                set(gca,'Color','k')  

                trial_offset = trial_offset + length(align_times);

                ylim([0 trial_offset+1])
                line([0,0],ylim,'linestyle','--','color',[.7 .7 .7]);
                line([500,500],ylim,'linestyle','--','color',[.55 .55 .55]);
                line([bins(1)*1000 bins(end)*1000],[trial_offset+1 trial_offset+1],'linestyle','--','color',tertile_color{contralateral_activity_tertile+cb});
                ylabel('sorted trial')
                
                max_fr = [max_fr max(psth_smooth)];
            end
            xlabel('Time from stim onset (ms)')
            subplot(length(stims_by_side{stim_side}),2,stim_num*2-1)
            xlabel('Time from stim onset (ms)')
            
            
            for stim_num = 1:length(stims_by_side{stim_side})
                subplot(length(stims_by_side{stim_side}),2,stim_num*2-1);
                ylim([-1 1.05*max(max_fr)])
                line([0,0],ylim,'linestyle','--','color',[.7 .7 1]);
                line([500 500],ylim,'linestyle','--','color',[.7 .7 1]);
            end
            
            if separate_contralateral_spontaneous && stim_num==1
                legend([h{1},h{2},h{3}],{'low','middle','high'},'TextColor','w')
            end
        end
        psth_fig.InvertHardcopy = 'off';
        saveas(psth_fig,[save_folder curr_figure_name '.tif'])
    end
    close all
end
end
