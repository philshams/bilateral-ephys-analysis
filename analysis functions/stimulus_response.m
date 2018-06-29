% -----------------------------------------------------------------
%                Stimulus Responses
% -----------------------------------------------------------------

%% PLOT PSTH BY STIMULUS AND OPTIONALLY BY CONTRALATERAL ACTIVITY TERTILE

% set up experiment parameters
tertile_color = {[.2 .2 .6],[.4 .4 .8],[.7 .7 1]};
contra_sites = [2 1];
avg_frs_ipsi_all = [];
avg_frs_all = [];
frs_ipsi_all = [];
frs_all = [];
h = {};
time_window = [-1,1.5];
psth_bin_size = 0.001;
stimulus_colors = [];
stim_color = [0 0 .3; .3 .3 .7; .4 .4 .9; .8 .8 .8; 1 1 1];
corr_fig = {}; hm = {};

% loop over probes
for site = sites

avg_frs_ipsi_by_side = [];
avg_frs_by_side = [];
frs_ipsi_by_side = [];
frs_side = [];
frs_ipsi_side_mean = [];
frs_side_mean = [];


% load data
spike_times_timeline = ephys_data.ephys_data{site,3};
cluster_IDs  = ephys_data.ephys_data{site,6};
num_neurons = length(cluster_IDs);

% load other probe's data if separating by contralateral activity
contra_site = contra_sites(site);
other_spike_times_timeline = ephys_data.ephys_data{contra_site,3};
other_cluster_IDs  = ephys_data.ephys_data{contra_site,6};
other_num_neurons = length(other_cluster_IDs);
if separate_contralateral_spontaneous    
    contra_text = [' by contra activity during ' num2str(epoch_to_split_by(1)*1000) '-' num2str(epoch_to_split_by(2)*1000) 'ms'];
else
    contra_text = '';
end

% look at contralateral stimulus presentation
stim_side = site;
    
    corr_fig{site} = figure('Position',[800 200 800 800]);  
    psth_fig = figure('Name',['site ' num2str(site) ' PSTH'],'Position',[700 50 800 1050]); hold on;
    stimulus_colors = [];
    
    for stim_num = 1:length(stims_by_side{stim_side})
        
            % align times by stimulus onset
            subplot(length(stims_by_side{stim_side}),1,stim_num); 
            stim = stims_by_side{stim_side}(stim_num);
            align_times_all = stimOn_times(ismember(stimIDs,stim));   
            stimulus_colors = [stimulus_colors; ones(length(align_times_all),3) .* stim_color(stim_num,:)];
            
            [frs, ~, ~, ~, ~,~,~] = ... 
            psth_and_smooth(spike_times_timeline, align_times_all, time_window, psth_bin_size, smooth_window, num_neurons);                
            [frs_ipsi, ~, ~, ~, ~,~,~] = ... 
            psth_and_smooth(other_spike_times_timeline, align_times_all, time_window, psth_bin_size, smooth_window, other_num_neurons); 
        
            [baseline_fr_ipsi, trial_tertiles] = get_avg_frs(other_spike_times_timeline, align_times_all, [-2 0], other_num_neurons);  
            [baseline_fr, ~] = get_avg_frs(spike_times_timeline, align_times_all, [-2 0], num_neurons);           
            
            baseline_fr_ipsi = mean(baseline_fr_ipsi);
            baseline_fr = mean(baseline_fr);

            % get trial nums of trial in current contralateral activity tertile      
            [avg_frs_ipsi, trial_tertiles_contra] = get_avg_frs(other_spike_times_timeline, align_times_all, epoch_to_split_by, other_num_neurons);  
            [avg_frs, trial_tertiles_ipsi] = get_avg_frs(spike_times_timeline, align_times_all, epoch_to_split_by, num_neurons);
    
            
        % loop across bins of magnitude of contralateral activity
        for contralateral_activity_tertile = 1:3^(separate_contralateral_spontaneous)
             
            if separate_contralateral_spontaneous
                align_times = align_times_all(trial_tertiles_contra==contralateral_activity_tertile);
                cb = 0;
            else
                align_times = align_times_all;
                cb = 2; % just adjusting the color of the PSTH
            end
            
            % get PSTH
            [psth_smooth_all, psth_smooth, psth_std_trials, bins, rasterX,rasterY,spikeCounts] = ... 
                psth_and_smooth(spike_times_timeline, align_times, time_window, psth_bin_size, smooth_window, num_neurons);

            % plot output
            fill([bins*1000; flip(bins*1000)], [(psth_smooth-psth_std_trials); flip(psth_smooth + psth_std_trials)], ...
                tertile_color{contralateral_activity_tertile+cb}, 'EdgeAlpha', 0, 'FaceAlpha', .3); 
            hold on
            h{contralateral_activity_tertile+cb} = plot(bins * 1000, psth_smooth, 'color',tertile_color{contralateral_activity_tertile+cb}, 'linewidth',2);

            title(['Population PSTH of ' site_sides{site} ' V1 to ' stim_sides{stim_side} ' ' contrasts{stim_num} ' stimulus' contra_text]);

            set(gca,'Color','k')
            ylim([0 35]);            
            line([0,0],ylim,'linestyle','--','color',[.7 .7 1]);
            line([500,500],ylim,'linestyle','--','color',[.5 .5 .8]);
            ylabel('avg firing rate');
            
            if separate_contralateral_spontaneous
                fill([epoch_to_split_by*1000 flip(epoch_to_split_by*1000)], [-100 -100 1000 1000], ...
                [.4 .4 .4], 'EdgeAlpha', 0, 'FaceAlpha', .07); 
            end
        end
        if separate_contralateral_spontaneous && stim_num==1
            legend([h{1},h{2},h{3}],{'low','middle','high'},'TextColor','w')
        end
        

        % concatenate data across the different stimuli
        avg_frs_ipsi_by_side = [avg_frs_ipsi_by_side; avg_frs_ipsi - baseline_fr_ipsi]; 
        avg_frs_by_side = [avg_frs_by_side; avg_frs - baseline_fr]; 

        figure(corr_fig{site}); hold on
        hm{stim_num} = scatter(mean(avg_frs_ipsi - baseline_fr_ipsi), mean(avg_frs - baseline_fr), 200, stim_color(stim_num,:), 'filled', 'MarkerEdgeColor', [.7 .7 .7]);
        figure(psth_fig)
        
        avg_frs_ipsi_all = [avg_frs_ipsi_all; zscore(avg_frs_ipsi)];
        avg_frs_all = [avg_frs_all; zscore(avg_frs)];

        frs_ipsi_by_side = [frs_ipsi_by_side; zscore(frs_ipsi)];
        frs_side = [frs_side; zscore(frs)];        
        
        frs_ipsi_all = [frs_ipsi_all; zscore(frs_ipsi)];
        frs_all = [frs_all; zscore(frs)];
            
        
    end
    xlabel('Time from stim onset (ms)')

% plot correlation for all trials on a given side
figure(corr_fig{site})
scatter(avg_frs_ipsi_by_side, avg_frs_by_side, 25, stimulus_colors, 'filled', 'MarkerEdgeColor', [.7 .7 .7])
set(gca,'Color','k')
xlabel('avg baseline-subtracted firing rate, ipsilateral to stimulus')
ylabel('avg baseline-subtracted firing rate, contralateral to stimulus')            
title(['population activity -- all ' stim_sides{site} ' stimuli,'...
         ' from ' num2str(epoch_to_split_by(1)*1000) '-' num2str(epoch_to_split_by(2)*1000) 'ms']);
     
% and plot best fit line
coefficients = polyfit(avg_frs_ipsi_by_side, avg_frs_by_side, 1);
xFit = linspace(min([avg_frs_ipsi_by_side; avg_frs_by_side])-3, max([avg_frs_ipsi_by_side; avg_frs_by_side])+3, 100);
% xFit = linspace(-15, 15, 10);
yFit = polyval(coefficients , xFit);

[R,P] = corrcoef(avg_frs_ipsi_by_side, avg_frs_by_side);

hold on;
plot(xFit, yFit, 'linestyle', '--', 'LineWidth', 2, 'color', [.5 .5 .5]);
plot(xFit, xFit, 'linestyle', ':', 'LineWidth', .5, 'color', [.9 .9 .9]);
xlim([-6 12])
ylim([-6 12])
text(10,-4,['r = ' num2str(round(R(2),2))],'color','white')
text(10,-5,['p = ' num2str(round(P(2),2))],'color','white')    

legend([hm{1},hm{2},hm{3},hm{4}],{'12%','25%','50%','100%'},'TextColor','w')




if plot_correlations

% plot correlation over time
[fr_correlations fr_correlations_P] = corrcoef([frs_side frs_ipsi_by_side]);
x_corr = fr_correlations(size(frs_all,2)+1:end, 1:size(frs_all,2));
auto_corr = fr_correlations(1:size(frs_all,2), 1:size(frs_all,2));
auto_corr_contra= fr_correlations(size(frs_all,2)+1:end, size(frs_all,2)+1:end);
x_corr_P = fr_correlations_P(size(frs_all,2)+1:end, 1:size(frs_all,2));
auto_corr_P = fr_correlations_P(1:size(frs_all,2), 1:size(frs_all,2));
auto_corr_contra_P = fr_correlations_P(size(frs_all,2)+1:end, size(frs_all,2)+1:end);

% ch = imagesc(fr_correlations);
figure('Position',[200 100 900 800]); hold on
ch = imagesc(auto_corr);
caxis([prctile(auto_corr(:),2) prctile(auto_corr(:),98)])
colormap(jet)
colorbar()
set(gca,'Color','k')

alpha_map = (1 - auto_corr_P).^50;
set(ch,'AlphaData',alpha_map);


title([site_sides{site} ' V1 -- auto correlation contralateral to stimulus (activity z-scored within condition)'])
set(ch, 'XData', bins*1000);
set(ch, 'YData', bins*1000);
plot([0 0],[bins(1)*1000 bins(end)*1000],'color',[.4 .4 .4],'linestyle','--')
plot([bins(1)*1000 bins(end)*1000],[0 0],'color',[.4 .4 .4],'linestyle','--')
plot([500 500],[bins(1)*1000 bins(end)*1000],'color',[.4 .4 .4],'linestyle','--')
plot([bins(1)*1000 bins(end)*1000],[500 500],'color',[.4 .4 .4],'linestyle','--')

xlim([bins(1)*1000 bins(end)*1000])
ylim([bins(1)*1000 bins(end)*1000])
xlabel('time (ms) contralateral to stimulus')
ylabel('time (ms) contralateral to stimulus')




figure('Position',[200 100 900 800]); hold on
ch = imagesc(auto_corr_contra);
caxis([prctile(auto_corr(:),2) prctile(auto_corr(:),98)])
colormap(jet)
colorbar()
set(gca,'Color','k')

alpha_map = (1 - auto_corr_contra_P).^100;
set(ch,'AlphaData',alpha_map);

title([site_sides{site} ' V1 -- auto correlation ipsilateral to stimulus (activity z-scored within condition)'])
set(ch, 'XData', bins*1000);
set(ch, 'YData', bins*1000);
plot([0 0],[bins(1)*1000 bins(end)*1000],'color',[.4 .4 .4],'linestyle','--')
plot([bins(1)*1000 bins(end)*1000],[0 0],'color',[.4 .4 .4],'linestyle','--')

plot([500 500],[bins(1)*1000 bins(end)*1000],'color',[.4 .4 .4],'linestyle','--')
plot([bins(1)*1000 bins(end)*1000],[500 500],'color',[.4 .4 .4],'linestyle','--')

xlim([bins(1)*1000 bins(end)*1000])
ylim([bins(1)*1000 bins(end)*1000])
xlabel('time (ms) ipsilateral to stimulus')
ylabel('time (ms) ipsilateral to stimulus')




figure('Position',[200 100 900 800]); hold on
ch = imagesc(x_corr);
% caxis([prctile(x_corr(:),5) prctile(x_corr(:),99)])
colormap(jet)
colorbar()
set(gca,'Color','k')

alpha_map = (1 - x_corr_P).^10;
set(ch,'AlphaData',alpha_map);


title([site_sides{site} ' V1 -- correlation of stimulus-evoked responses and contralateral activity (activity z-scored within condition)'])
set(ch, 'XData', bins*1000);
set(ch, 'YData', bins*1000);
plot([0 0],[bins(1)*1000 bins(end)*1000],'color',[.4 .4 .4],'linestyle','--')
plot([bins(1)*1000 bins(end)*1000],[0 0],'color',[.4 .4 .4],'linestyle','--')
plot([500 500],[bins(1)*1000 bins(end)*1000],'color',[.4 .4 .4],'linestyle','--')
plot([bins(1)*1000 bins(end)*1000],[500 500],'color',[.4 .4 .4],'linestyle','--')
xlim([bins(1)*1000 bins(end)*1000])
ylim([bins(1)*1000 bins(end)*1000])
xlabel('time (ms) contralateral to stimulus')
ylabel('time (ms) ipsilateral to stimulus')

end



end




% plot correlation for all trials
figure('Position',[800 200 800 800]);  
scatter(avg_frs_ipsi_all, avg_frs_all, 15, 'filled' ,'MarkerFaceColor',tertile_color{2}, 'MarkerEdgeColor', [.7 .7 .7])
set(gca,'Color','k')
xlabel('zscore over trials, ipsilateral to stimulus')
ylabel('zscore over trials, contralateral to stimulus')            
title(['population activity -- all trials,' ...
         ' from ' num2str(epoch_to_split_by(1)*1000) '-' num2str(epoch_to_split_by(2)*1000) 'ms']);

% and plot best fit line
coefficients = polyfit(avg_frs_ipsi_all, avg_frs_all, 1);
xFit = linspace(-4, 4, 100);
yFit = polyval(coefficients , xFit);

[R,P] = corrcoef(avg_frs_ipsi_all, avg_frs_all);

hold on;
plot(xFit, yFit, 'linestyle', '--', 'LineWidth', 2, 'color', [.5 .5 .5]);
plot(xFit, xFit, 'linestyle', ':', 'LineWidth', .5, 'color', [.9 .9 .9]);
xlim([-4 4])
ylim([-4 4])
text(3,-2.5,['r = ' num2str(round(R(2),2))],'color','white')
text(3,-3,['p = ' num2str(round(P(2),2))],'color','white')   
