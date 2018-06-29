% -----------------------------------------------------------------
%                           Depth Plot of Correlations
% -----------------------------------------------------------------



if strcmp(use_stim_type,'weak') || strcmp(use_stim_type,'low')
    curr_baseline_correlations_stim = baseline_correlations_stim_weak;
    stim_type = ' weak';
elseif strcmp(use_stim_type,'strong') || strcmp(use_stim_type,'high')
    curr_baseline_correlations_stim = baseline_correlations_stim_strong;
    stim_type = ' strong';
else
    curr_baseline_correlations_stim = baseline_correlations_stim;
    stim_type = ' all';
end

stim_words = [stim_type ' stimulus - baseline'];



curr_baseline_correlations_stim{1,1}(:,4) = 0;
curr_baseline_correlations_stim{1,2}(:,4) = 0;    


% baseline_correlations{1,2}(:,4) % p value
% baseline_correlations{1,2}(:,5) % mean baseline fr

% load right data
right_spike_times_timeline = ephys_data.ephys_data{1,3};
right_spike_clusters = ephys_data.ephys_data{1,5};
right_cluster_IDs  = ephys_data.ephys_data{1,6};
right_templateDepths = ephys_data.ephys_data{1,8};
right_templateDepths(6) = 980;

% load left data
left_spike_times_timeline = ephys_data.ephys_data{2,3};
left_spike_clusters = ephys_data.ephys_data{2,5};
left_cluster_IDs  = ephys_data.ephys_data{2,6};
left_templateDepths = ephys_data.ephys_data{2,8};

% mean / std
[right_depths, right_depths_ind] = sort(right_templateDepths);
right_corr_sorted = (curr_baseline_correlations_stim{1,1}(right_depths_ind,2) - baseline_correlations{1,1}(right_depths_ind,2));
smooth_depths_right = movmean(right_depths,8);
smooth_corr_right = movmean(right_corr_sorted,8);

[left_depths, left_depths_ind] = sort(left_templateDepths);
left_corr_sorted = (curr_baseline_correlations_stim{1,2}(left_depths_ind,2) - baseline_correlations{1,2}(left_depths_ind,2));
smooth_depths_left = movmean(left_depths,8);
smooth_corr_left = movmean(left_corr_sorted,8);


    % r < .05
right_modulated = ones(length(curr_baseline_correlations_stim{1,1}(:,2)),1);
    % r > .05
right_modulated( (curr_baseline_correlations_stim{1,1}(:,2) - baseline_correlations{1,1}(:,2)) > 0 ) = 2;
    % r > .1
right_modulated( (curr_baseline_correlations_stim{1,1}(:,2) - baseline_correlations{1,1}(:,2)) > .03 ) = 3;
    % r > .2
right_modulated( (curr_baseline_correlations_stim{1,1}(:,2) - baseline_correlations{1,1}(:,2)) > .1 ) = 4;
    % stim unresponsive
right_modulated( curr_baseline_correlations_stim{1,1}(:,4) > .05 ) = right_modulated( curr_baseline_correlations_stim{1,1}(:,4) > .05 ) + 4;
    % by fr

right_modulated( (baseline_correlations{1,1}(:,5) > 2) & (baseline_correlations{1,1}(:,5) < 10) ) = ...
    right_modulated( (baseline_correlations{1,1}(:,5) > 2) & (baseline_correlations{1,1}(:,5) < 10) ) + 8;

right_modulated( (baseline_correlations{1,1}(:,5) >= 10) ) = ...
    right_modulated( (baseline_correlations{1,1}(:,5) >= 10) ) + 16;


    % r < .05
left_modulated = ones(length(curr_baseline_correlations_stim{1,2}(:,2)),1);
    % r > .05
left_modulated( (curr_baseline_correlations_stim{1,2}(:,2) - baseline_correlations{1,2}(:,2)) > 0 ) = 2;
    % r > .1
left_modulated( (curr_baseline_correlations_stim{1,2}(:,2) - baseline_correlations{1,2}(:,2)) > .03 ) = 3;
    % r > .2
left_modulated( (curr_baseline_correlations_stim{1,2}(:,2) - baseline_correlations{1,2}(:,2)) > .1 ) = 4;
    % stim unresponsive
left_modulated( curr_baseline_correlations_stim{1,2}(:,4) > .05 ) = left_modulated( curr_baseline_correlations_stim{1,2}(:,4) > .05 ) + 4;
    % by fr
if show_baseline_fr
    left_modulated( (baseline_correlations{1,2}(:,5) > 2) & (baseline_correlations{1,2}(:,5) < 10) ) = ...
        left_modulated( (baseline_correlations{1,2}(:,5) > 2) & (baseline_correlations{1,2}(:,5) < 10) ) + 8;

    left_modulated( (baseline_correlations{1,2}(:,5) >= 10) ) = ...
        left_modulated( (baseline_correlations{1,2}(:,5) >= 10) ) + 16;
end



celltype_labels = {'         \Deltar < 0','0.0 < \Deltar < .03','.03 < \Deltar < 0.1', '0.1 < \Deltar', '',...
         'low baseline fr (< 2 Hz)', 'med baseline fr (> 2 Hz)', 'high baseline fr (> 10 Hz)'};


category_colors = repmat({[.4 .4 .6 ], [.5 .5 .5], [.8 .8 .5 ], [1 1 0 ]},1,6) ;
category_sizes = {6,6,6,6,6,6,6,6,10,10,10,10,10,10,10,10,14,14,14,14,14,14,14,14};             
category_markers = repmat({'o','o','o','o','s','s','s','s'},1,3);

f = figure('Position',[2003 -316 1000 1307],'Name',['Unit Depths Stim ' num2str(stim_side)]); hold on
opt.xyOri = 'normal';

% plot left
subplot(1,2,1); set(gca,'Color','k'); ylim([-200, 1100]); set(gca,'YDir','reverse'); 
use_colors = unique(left_modulated);

plotSpread_correlations(left_templateDepths,'categoryIdx', left_modulated,...
                'categoryColors',category_colors(use_colors),'categoryMarkers',category_markers(use_colors),'categorySizes',category_sizes(use_colors));

ylabel('Depth (\mum)'); title(['left V1 correlations to contra' stim_words])
set(gca,'XTick',[])
x_limit = xlim;
xlim([-.2 x_limit(2)])
plot(smooth_corr_left*5, smooth_depths_left, 'linewidth',2,'color','r');

% plot right
subplot(1,2,2); set(gca,'Color','k'); ylim([-200, 1100])
set(gca,'YDir','reverse'); title(['right V1 correlations to contra' stim_words])
use_colors = unique(right_modulated);

plotSpread_correlations(right_templateDepths,'categoryIdx', right_modulated,...
                'categoryColors',category_colors(use_colors),'categoryMarkers',category_markers(use_colors),'categorySizes',category_sizes(use_colors));

set(gca,'XTick',[]);
xlim([-.2 x_limit(2)])
plot(smooth_corr_right*5, smooth_depths_right, 'linewidth',2,'color','r');


% add legend
p1 = plot(-10,-10,'marker',category_markers{1},'markers', category_sizes{10},'color',category_colors{1},...
                        'lineStyle','none', 'linewidth',1,'MarkerFaceColor',category_colors{1});
p2 = plot(-10,-10,'marker',category_markers{1},'markers', category_sizes{10},'color',category_colors{2},...
                        'lineStyle','none', 'linewidth',1,'MarkerFaceColor',category_colors{2});
p3 = plot(-10,-10,'marker',category_markers{1},'markers', category_sizes{10},'color',category_colors{3},...
                        'lineStyle','none', 'linewidth',1,'MarkerFaceColor',category_colors{3});
p4 = plot(-10,-10,'marker',category_markers{1},'markers', category_sizes{10},'color',category_colors{4},...
                        'lineStyle','none', 'linewidth',1,'MarkerFaceColor',category_colors{4});
                   
p7 = plot(-10,-10,'marker',category_markers{1},'markers', category_sizes{1},'color',category_colors{3},...
                        'lineStyle','none', 'linewidth',1,'MarkerFaceColor',category_colors{3});
p8 = plot(-10,-10,'marker',category_markers{1},'markers', category_sizes{10},'color',category_colors{3},...
                        'lineStyle','none', 'linewidth',1,'MarkerFaceColor',category_colors{3});
p9 = plot(-10,-10,'marker',category_markers{1},'markers', category_sizes{20},'color',category_colors{3},...
                        'lineStyle','none', 'linewidth',1,'MarkerFaceColor',category_colors{3});
pblack1 = plot(-10,-10,'marker',category_markers{1},'markers', category_sizes{1},'color',[0 0 0],...
                        'lineStyle','none', 'linewidth',1,'MarkerFaceColor',[0 0 0]);
pblack2 = plot(-10,-10,'marker',category_markers{1},'markers', category_sizes{1},'color',[0 0 0],...
                        'lineStyle','none', 'linewidth',1,'MarkerFaceColor',[0 0 0]);                    
                    

legend([p1 p2 p3 p4 pblack2 p7 p8 p9 ],celltype_labels,'location','north','TextColor','w');



