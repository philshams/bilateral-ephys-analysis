% -----------------------------------------------------------------
%          Depth Plot of Stim-Responsive / Contra-Modulated Units
% -----------------------------------------------------------------



% load right data
right_spike_times_timeline = ephys_data.ephys_data{1,3};
right_spike_clusters = ephys_data.ephys_data{1,5};
right_cluster_IDs  = ephys_data.ephys_data{1,6};
right_templateDepths = ephys_data.ephys_data{1,8};

% load left data
left_spike_times_timeline = ephys_data.ephys_data{2,3};
left_spike_clusters = ephys_data.ephys_data{2,5};
left_cluster_IDs  = ephys_data.ephys_data{2,6};
left_templateDepths = ephys_data.ephys_data{2,8};

% loop left then right stims
for stim_side = 1:2

    % get stim times
    stims = stims_by_side{stim_side};
    align_times_all = stimOn_times(ismember(stimIDs,stims));      

    % assign units as stim-responsive or not (get p-values) // assign units as modulated by population activity (option for contra or ipsi)
        % loop through to z-score for each stimulus
        pop_fr_right = [];
        pop_fr_left = [];
    for stim_num = 1:length(stims_by_side{stim_side})
        
        stim = stims_by_side{stim_side}(stim_num);
        align_times_stim = stimOn_times(ismember(stimIDs,stim));          
        
        [pop_stim_fr_right, ~] = get_avg_frs(right_spike_times_timeline, align_times_stim, [0 .5], length(right_cluster_IDs));  
        [pop_stim_fr_left, ~] = get_avg_frs(left_spike_times_timeline, align_times_stim, [0 .5], length(left_cluster_IDs));     
        
        pop_fr_right = [pop_fr_right; zscore(pop_stim_fr_right)];
        pop_fr_left = [pop_fr_left; zscore(pop_stim_fr_left)];  
        
    end
    
    right_p_value = zeros(length(right_cluster_IDs),1);
    left_p_value = zeros(length(left_cluster_IDs),1);    

    right_pop_corr = zeros(length(right_cluster_IDs),2);
    left_pop_corr = zeros(length(left_cluster_IDs),2);       
    
    right_baseline_fr = zeros(length(right_cluster_IDs),1);
    left_baseline_fr = zeros(length(left_cluster_IDs),1);        
    
    for unit_num = 1:length(right_cluster_IDs)
        unit = right_cluster_IDs(unit_num);
        [baseline_fr_right, ~] = get_avg_frs(right_spike_times_timeline(right_spike_clusters==unit), align_times_all, [-1 0], 1);  
        [stim_fr_right, ~] = get_avg_frs(right_spike_times_timeline(right_spike_clusters==unit), align_times_all, [0 .5], 1);          
        [h, p] = ttest(stim_fr_right, baseline_fr_right, 'Tail', 'Right');
        right_p_value(unit_num) = p;
        right_baseline_fr(unit_num) = mean(baseline_fr_right);
        
        % loop through to z-score for each stimulus
        unit_fr_right = [];
        for stim_num = 1:length(stims_by_side{stim_side})
            stim = stims_by_side{stim_side}(stim_num);
            align_times_stim = stimOn_times(ismember(stimIDs,stim));       
            [stim_fr_right, ~] = get_avg_frs(right_spike_times_timeline(right_spike_clusters==unit), align_times_stim, [0 .5], 1);  
            unit_fr_right = [unit_fr_right; zscore(stim_fr_right)];
            
        end
        
        [R, P] = corrcoef(pop_fr_right,unit_fr_right);
        right_pop_corr(unit_num,1) = R(2);
        right_pop_corr(unit_num,2) = P(2);
        
    end
    
    for unit_num = 1:length(left_cluster_IDs)
        unit = left_cluster_IDs(unit_num);
        [baseline_fr_left, ~] = get_avg_frs(left_spike_times_timeline(left_spike_clusters==unit), align_times_all, [-1 0], 1);        
        [stim_fr_left, ~] = get_avg_frs(left_spike_times_timeline(left_spike_clusters==unit), align_times_all, [0 .5], 1);               
        [h, p] = ttest(stim_fr_left, baseline_fr_left, 'Tail', 'Right');
        left_p_value(unit_num) = p;
        left_baseline_fr(unit_num) = mean(baseline_fr_left);
        
        % loop through to z-score for each stimulus
        unit_fr_left = [];
        for stim_num = 1:length(stims_by_side{stim_side})

            stim = stims_by_side{stim_side}(stim_num);
            align_times_stim = stimOn_times(ismember(stimIDs,stim));  
            [stim_fr_left, ~] = get_avg_frs(left_spike_times_timeline(left_spike_clusters==unit), align_times_stim, [0 .5], 1);  
            unit_fr_left = [unit_fr_left; zscore(stim_fr_left)];        
            
        end        
        
        [R, P] = corrcoef(pop_fr_left,unit_fr_left);
        left_pop_corr(unit_num,1) = R(2);
        left_pop_corr(unit_num,2) = P(2);        
        
    end
   
   if stim_side == 1
        baseline_correlations_stim{1,1}(:,4) = right_p_value;
        baseline_correlations{1,1}(:,4) = right_p_value;
        baseline_correlations_stim{1,1}(:,5) = right_baseline_fr;
        baseline_correlations{1,1}(:,5) = right_baseline_fr;        
   elseif stim_side == 2
        baseline_correlations_stim{1,2}(:,4) = left_p_value;
        baseline_correlations{1,2}(:,4) = left_p_value;
        baseline_correlations_stim{1,2}(:,5) = left_baseline_fr;
        baseline_correlations{1,2}(:,5) = left_baseline_fr;        
   end
   
   
    % Plot cell type by depth
        % non responsive, non-modulated
    left_modulated = ones(length(left_pop_corr),1);
        % non responsive, contra-modulated
    left_modulated( (left_p_value >= .05) & (left_pop_corr(:,2) < .05) & (left_pop_corr(:,1) > .1) ) = 2;
        % responsive, non-modulated
    left_modulated( (left_p_value < .05) ) = 3;
        % responsive, contra-modulated        
    left_modulated( (left_p_value < .05) & (left_pop_corr(:,2) < .05) & (left_pop_corr(:,1) > .1) ) = 4;
    
        % non responsive, non-modulated
    right_modulated = ones(length(right_pop_corr),1);
        % non responsive, contra-modulated
    right_modulated( (right_p_value >= .05) & (right_pop_corr(:,2) < .05) & (right_pop_corr(:,1) > .1) ) = 2;
        % responsive, non-modulated
    right_modulated( (right_p_value < .05) ) = 3;
        % responsive, contra-modulated        
    right_modulated( (right_p_value < .05) & (right_pop_corr(:,2) < .05) & (right_pop_corr(:,1) > .1) ) = 4;
    
        
%     
%     
%     celltype_labels = {'stim-unresponsive / independent of ipsi activity','stim-unresponsive / modulated by ipsi activity'...
%                                 'stim-responsive / independent of ipsi activity','stim-responsive / modulated by ipsi activity'};
%     
%     category_colors = {[1 1 1], [1 1 1],[.4 .4 1],[.4 .4 1],};
%     category_markers = {'.','*','.','*'};
%     
%     f = figure('Position',[727 42 1039 1074],'Name',['Unit Depths Stim ' num2str(stim_side)]); hold on
%     opt.xyOri = 'normal';
%     
%     % plot left
%     subplot(1,2,1); set(gca,'Color','k'); ylim([-200, 1100]); set(gca,'YDir','reverse'); 
%     use_colors = unique(left_modulated);
%     
%     plotSpread(left_templateDepths,'categoryIdx', left_modulated,...
%                     'categoryColors',category_colors(use_colors),'categoryMarkers',category_markers(use_colors));
%                 
%     ylabel('Depth (\mum)'); title(['left hemisphere / ' stim_sides{stim_side} ' visual stimulus'])
%     set(gca,'XTick',[])
%     
%     % plot right
%     subplot(1,2,2); set(gca,'Color','k'); ylim([-200, 1100])
%     set(gca,'YDir','reverse'); title(['right hemisphere / ' stim_sides{stim_side} ' visual stimulus'])
%     use_colors = unique(right_modulated);
%     
%     plotSpread(right_templateDepths,'categoryIdx', right_modulated,...
%                     'categoryColors',category_colors(use_colors),'categoryMarkers',category_markers(use_colors));
%                 
%     set(gca,'XTick',[]);
%     
%     % add legend
%     subplot(1,2,(stim_side-2)*-1+1)
%     legend(celltype_labels,'location','north','TextColor','w');
%     

end