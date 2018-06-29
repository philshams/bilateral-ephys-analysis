
% -----------------------------------------------------------------
%               Pairwise Correlations
% -----------------------------------------------------------------



corr_window = 0.1; % duration to measure correlation
epoch_start_times = -1.25:.025:1.65;

show_correlation_plot = false

mean_R_contra_array = zeros(2,length(epoch_start_times));
std_R_contra_array = zeros(2,length(epoch_start_times));
frac_sig_R_contra_array = zeros(2,length(epoch_start_times));
mean_R_ipsi_array = zeros(4,length(epoch_start_times));
std_R_ipsi_array = zeros(4,length(epoch_start_times));
frac_sig_R_ipsi_array = zeros(4,length(epoch_start_times));


for epoch_num = 1:length(epoch_start_times)

    epoch_start = epoch_start_times(epoch_num)
    epoch_to_correlate = [epoch_start epoch_start+corr_window];    

    mean_R_contra = zeros(2,1);
    std_R_contra = zeros(2,1);
    num_sig_R_contra = zeros(2,1);
    mean_R_ipsi = zeros(4,1);
    std_R_ipsi = zeros(4,1);
    num_sig_R_ipsi = zeros(4,1);
    
    unit_correlation
    mean_R_ipsi_array(:,epoch_num) = mean_R_ipsi;
    std_R_ipsi_array(:,epoch_num) = std_R_ipsi;
    frac_sig_R_ipsi_array(:,epoch_num) = frac_sig_R_ipsi;
    
    contralateral_unit_correlation
    mean_R_contra_array(:,epoch_num) = mean_R_contra;
    std_R_contra_array(:,epoch_num) = std_R_contra;
    frac_sig_R_contra_array(:,epoch_num) = frac_sig_R_contra;

end

close all

for site = 1:2

    stim_colors = {[.75 0 .75],[0 .7 0]}; %left right
    site_ind = site + (site-1);
    
    % plot right probe ipsi corr
    f = figure('Position',[669 339 1068 574]); hold on
    set(gca,'Color','k') 

    time_axis = (epoch_start_times + corr_window / 2)*1000;

    p1 = plot(time_axis, mean_R_ipsi_array(site_ind,:),'color',stim_colors{1},'linewidth',2);
    fill([time_axis flip(time_axis)], [(mean_R_ipsi_array(site_ind,:)-std_R_ipsi_array(site_ind,:)) flip(mean_R_ipsi_array(site_ind,:)+std_R_ipsi_array(site_ind,:))], ...
        stim_colors{1}, 'EdgeAlpha', 0, 'FaceAlpha', .3); 

    p2 = plot(time_axis, mean_R_ipsi_array(site_ind+1,:),'color',stim_colors{2},'linewidth',2);
    fill([time_axis flip(time_axis)], [(mean_R_ipsi_array(site_ind+1,:)-std_R_ipsi_array(site_ind+1,:)) flip(mean_R_ipsi_array(site_ind+1,:)+std_R_ipsi_array(site_ind+1,:))], ...
        stim_colors{2}, 'EdgeAlpha', 0, 'FaceAlpha', .3); 

    ylim(ylim)

    line([0,0],ylim,'linestyle','--','color',[.7 .7 .7]);
    line([500,500],ylim,'linestyle','--','color',[.55 .55 .55]);

    title([site_sides{site} ' V1 avg ipsi pairwise correlation coefficients over time'])
    legend([p1 p2],{'left stim','right stim'},'TextColor','w')
    ylabel('correlation coefficient')
    xlabel('time (ms)')

    
    
    % plot right probe ipsi num sig corr
    f = figure('Position',[669 339 1068 574]); hold on
    set(gca,'Color','k') 

    p1 = plot(time_axis, frac_sig_R_ipsi_array(site_ind,:),'color',stim_colors{1},'linewidth',3);
    p2 = plot(time_axis, frac_sig_R_ipsi_array(site_ind+1,:),'color',stim_colors{2},'linewidth',3);

    ylim(ylim)
    line([0,0],ylim,'linestyle','--','color',[.7 .7 .7]);
    line([500,500],ylim,'linestyle','--','color',[.55 .55 .55]);

    title([site_sides{site} ' V1 fraction significant ipsi pairwise correlations over time'])
    legend([p1 p2],{'left stim','right stim'},'TextColor','w')
    ylabel('franction significant pairwise noise correlations')
    xlabel('time (ms)')

end




    stim_colors = {[.75 0 .75],[0 .7 0]}; %left right
    
    % plot right probe ipsi corr
    f = figure('Position',[669 339 1068 574]); hold on
    set(gca,'Color','k') 

    time_axis = (epoch_start_times + corr_window / 2)*1000;

    p1 = plot(time_axis, mean_R_contra_array(1,:),'color',stim_colors{1},'linewidth',2);
    fill([time_axis flip(time_axis)], [(mean_R_contra_array(1,:)-std_R_contra_array(1,:)) flip(mean_R_contra_array(1,:)+std_R_contra_array(1,:))], ...
        stim_colors{1}, 'EdgeAlpha', 0, 'FaceAlpha', .3); 

    p2 = plot(time_axis, mean_R_contra_array(2,:),'color',stim_colors{2},'linewidth',2);
    fill([time_axis flip(time_axis)], [(mean_R_contra_array(2,:)-std_R_contra_array(2,:)) flip(mean_R_contra_array(2,:)+std_R_contra_array(2,:))], ...
        stim_colors{2}, 'EdgeAlpha', 0, 'FaceAlpha', .3); 

    ylim(ylim)

    line([0,0],ylim,'linestyle','--','color',[.7 .7 .7]);
    line([500,500],ylim,'linestyle','--','color',[.55 .55 .55]);

    title(['avg pairwise contralateral correlation coefficients over time'])
    legend([p1 p2],{'left stim','right stim'},'TextColor','w')
    ylabel('correlation coefficient')
    xlabel('time (ms)')
    
    
    
    % plot right probe contra num sig corr
    f = figure('Position',[669 339 1068 574]); hold on
    set(gca,'Color','k') 

    p1 = plot(time_axis, frac_sig_R_contra_array(1,:),'color',stim_colors{1},'linewidth',3);
    p2 = plot(time_axis, frac_sig_R_contra_array(2,:),'color',stim_colors{2},'linewidth',3);

    ylim(ylim)
    line([0,0],ylim,'linestyle','--','color',[.7 .7 .7]);
    line([500,500],ylim,'linestyle','--','color',[.55 .55 .55]);

    title(['fraction significant contralateral pairwise correlations over time'])
    legend([p1 p2],{'left stim','right stim'},'TextColor','w')
    ylabel('franction significant pairwise noise correlations')
    xlabel('time (ms)')
    xlim([-1250 1750])
    
    
    