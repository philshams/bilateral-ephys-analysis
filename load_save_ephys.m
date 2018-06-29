%% Load ephys data (single long recording)

ephys_data = cell(length(sites),7);

for site = sites
    
[ephys_path,ephys_exists] = AP_cortexlab_filename(animal,day,experiment,'ephys',site);

if ephys_exists && load_parts.ephys
    
    if verbose; disp('Loading ephys...'); end;
    
    acqLive_channel = 2;
    load_lfp = true;
    
    % Load clusters, if they exist
    cluster_filename = [ephys_path filesep 'cluster_groups.csv'];
    if exist(cluster_filename,'file')
        fid = fopen(cluster_filename);
        cluster_groups = textscan(fid,'%d%s','HeaderLines',1);
        fclose(fid);
    end
    
    % Apparently now sometimes it's a different filename/type, if that
    % exists overwrite the other one
    cluster_filename = [ephys_path filesep 'cluster_group.tsv'];
    if exist(cluster_filename,'file')
        fid = fopen(cluster_filename);
        cluster_groups = textscan(fid,'%d%s','HeaderLines',1);
        fclose(fid);
    end
    
    % Load sync/photodiode
    load(([ephys_path filesep 'sync.mat']));
    
    % Read header information
    header_path = [ephys_path filesep 'dat_params.txt'];
    header_fid = fopen(header_path);
    header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
    fclose(header_fid);
    
    header = struct;
    for i = 1:length(header_info{1})
        header.(header_info{1}{i}) = header_info{2}{i};
    end
    
    % Load spike data
    if isfield(header,'sample_rate')
        ephys_sample_rate = str2num(header.sample_rate);
    elseif isfield(header,'ap_sample_rate')
        ephys_sample_rate = str2num(header.ap_sample_rate);
    end
    spike_times = double(readNPY([ephys_path filesep 'spike_times.npy']))./ephys_sample_rate;
    spike_templates = readNPY([ephys_path filesep 'spike_templates.npy']);
    spike_clusters = uint32(readNPY([ephys_path filesep 'spike_clusters.npy']));
    templates = readNPY([ephys_path filesep 'templates.npy']);
    channel_positions = readNPY([ephys_path filesep 'channel_positions.npy']);
    channel_map = readNPY([ephys_path filesep 'channel_map.npy']);
    winv = readNPY([ephys_path filesep 'whitening_mat_inv.npy']);
    template_amplitudes = readNPY([ephys_path filesep 'amplitudes.npy']);
    
    % Flip channel map and positions if banks are reversed
    % (this was only for phase 2, so setting false by default)
    flipped_banks = false;
    if flipped_banks
        channel_map = [channel_map(61:end);channel_map(1:60)];
        channel_positions = [channel_positions(61:end,:);channel_positions(1:60,:)];
    end
    
    % Default channel map/positions are from end: make from surface
    channel_positions(:,2) = max(channel_positions(:,2)) - channel_positions(:,2);
    
    % Load LFP
    n_channels = str2num(header.n_channels);
    %lfp_filename = [ephys_path filesep 'lfp.dat']; (this is old)
    [data_path,data_path_exists] = AP_cortexlab_filename(animal,day,experiment,'ephysraw',site);
    lfp_dir = dir([data_path 'experiment*-1_0.dat']);
    lfp_filename = [data_path lfp_dir.name];
    if load_lfp && exist(lfp_filename,'file')
        lfp_sample_rate = str2num(header.lfp_sample_rate);
        lfp_cutoff = str2num(header.filter_cutoff);
        
        fid = fopen(lfp_filename);
        % define where/how much of LFP to load
        lfp_skip_minutes = 10; % move to N minutes after recording start
        lfp_load_start = (lfp_sample_rate*60*lfp_skip_minutes*n_channels);
        lfp_load_samples = 1e6;
        % load LFP
        fseek(fid,lfp_load_start,'bof');
        lfp_all = fread(fid,[n_channels,lfp_load_samples],'int16'); % pull snippet
        fclose(fid);
        % eliminate non-connected channels
        lfp = lfp_all(channel_map+1,:);
        clear lfp_all;
        
        lfp_t = [(lfp_load_start/n_channels):(lfp_load_start/n_channels)+lfp_load_samples-1]/lfp_sample_rate;
    end
    
    % Get acqLive times for current experiment
    experiment_ephys_starts = sync(acqLive_channel).timestamps(sync(acqLive_channel).values == 1);
    experiment_ephys_stops = sync(acqLive_channel).timestamps(sync(acqLive_channel).values == 0);
    
    % (get folders with only a number - those're the experiment folders)
    experiments_dir = dir(AP_cortexlab_filename(animal,day,experiment,'expInfo'));
    experiments_num_idx = cellfun(@(x) ~isempty(x), regexp({experiments_dir.name},'^\d*$'));
    experiment_num = experiment == cellfun(@str2num,{experiments_dir(experiments_num_idx).name});
    acqlive_ephys_currexpt = [experiment_ephys_starts(experiment_num), ...
        experiment_ephys_stops(experiment_num)];
    
    % Get the spike/lfp times in timeline time (accounts for clock drifts)
    spike_times_timeline = AP_clock_fix(spike_times,acqlive_ephys_currexpt,acqLive_timeline);
    if load_lfp && exist(lfp_filename,'file')
        lfp_t_timeline = AP_clock_fix(lfp_t,acqlive_ephys_currexpt,acqLive_timeline);
    end
    
    % Get the depths of each template
    % (by COM - this used to not work but now looks ok)
    [spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
        templatePositionsAmplitudes(templates,winv,channel_positions(:,2),spike_templates,template_amplitudes);
    
    % Get the waveform duration of all templates (channel with largest negative amp)
%     [~,max_site] = max(max(abs(templates),[],2),[],3);
     [~,max_site] = min(min(templates,[],2),[],3);
    templates_max = nan(size(templates,1),size(templates,2));
    for curr_template = 1:size(templates,1)
        templates_max(curr_template,:) = ...
            templates(curr_template,:,max_site(curr_template));
    end
    waveforms = templates_max;
    
    % Get trough-to-peak time for each template
    templates_max_signfix = bsxfun(@times,templates_max, ...
        sign(abs(min(templates_max,[],2)) - abs(max(templates_max,[],2))));
    
    [~,waveform_trough] = min(templates_max,[],2);
    [~,waveform_peak_rel] = arrayfun(@(x) ...
        max(templates_max(x,waveform_trough(x):end),[],2), ...
        transpose(1:size(templates_max,1)));
    waveform_peak = waveform_peak_rel + waveform_trough;
    
    templateDuration = waveform_peak - waveform_trough;
    templateDuration_us = (templateDuration/ephys_sample_rate)*1e6;
    
    % Eliminate spikes that were classified as not "good"
    if exist('cluster_groups','var')
        
        if verbose; 
            if include_MUA; disp('Cleansing of noise...');
            else; disp('Cleansing of noise and MUA...'); end; end;
        
        if include_MUA
        good_templates_idx = uint32(cluster_groups{1}( ...
            strcmp(cluster_groups{2},'good') | strcmp(cluster_groups{2},'mua')));
        else
        good_templates_idx = uint32(cluster_groups{1}( ...
            strcmp(cluster_groups{2},'good')));            
        end
        
        % for merged and split units, use the higest contributing template
        % as the representative 'good' template
        new_clusters =  good_templates_idx(good_templates_idx > size(templates,1));
        template_clusters =  good_templates_idx(good_templates_idx <= size(templates,1));
        good_cluster_templates_idx = zeros(size(good_templates_idx), 'uint32');
        good_cluster_templates_idx(1:length(template_clusters)) = template_clusters;
        
        for cluster = 1:length(new_clusters)
            [N, edges] = histcounts(spike_templates(spike_clusters==new_clusters(cluster)), 0:2000);
            [~, max_contrib_ind] = max(N);
            good_cluster_templates_idx(length(template_clusters) + cluster) = max_contrib_ind-1;
        end
        good_cluster_templates_idx = good_cluster_templates_idx + 1;
        
        % Throw out all non-good template data
        if max_depth_to_analyze
            templateDepths_in_order = templateDepths(good_cluster_templates_idx);
            templateDepths_threshold = templateDepths_in_order < (max_depth_to_analyze + (3840-insertion_depth));        
            good_cluster_templates_idx = good_cluster_templates_idx(templateDepths_threshold);
            good_templates_idx = good_templates_idx(templateDepths_threshold);
        end
        
        templates = templates(good_cluster_templates_idx,:,:);
        templateDepths = templateDepths(good_cluster_templates_idx);
        waveforms = waveforms(good_cluster_templates_idx,:);
        templateDuration = templateDuration(good_cluster_templates_idx);
        templateDuration_us = templateDuration_us(good_cluster_templates_idx);
        
        
        % Throw out all non-good spike data
        good_spike_idx = ismember(spike_clusters,good_templates_idx);
        spike_times = spike_times(good_spike_idx);
        spike_clusters = spike_clusters(good_spike_idx);
        template_amplitudes = template_amplitudes(good_spike_idx);
        spikeDepths = spikeDepths(good_spike_idx);
        spike_times_timeline = spike_times_timeline(good_spike_idx);
        
        
    elseif ~exist('cluster_groups','var')
        if verbose; disp('Clusters not yet sorted'); end;
    end
    
end


%% Classify spikes

if ephys_exists && load_parts.ephys && exist('cluster_groups','var')
    if verbose; disp('Classifying spikes...'); end;
    

      
    % Cortical classification (like Bartho JNeurophys 2004)
    waveform_duration_cutoff = 400;
    narrow = templateDuration_us <= waveform_duration_cutoff;
    wide = templateDuration_us > waveform_duration_cutoff;
        
    waveform_t = 1e3*((0:size(templates,2)-1)/ephys_sample_rate);
    
    if verbose
        
        % Plot the waveforms and spike statistics
        figure('Position',[600 400 600 400],'Name',['Waveforms Site ' num2str(site)]);
        
        subplot(1,1,1); hold on;
        p = plot(waveform_t,waveforms');
        set(gca,'Color','k')
        set(p(wide),'color','w')
        set(p(narrow),'color','r')
        xlabel('Time (ms)')
        title('Cortical Waveforms');
        xlim([.75 2.75]);
        legend([p(find(wide,1)),p(find(narrow,1))],{'\color{white} Wide','\color{red} Narrow'})

        
        
        % Plot cell type by depth
        celltype_labels = {'Wide Spike','Narrow Spike'};
        celltypes = wide.*1 + narrow.*2 ;
        use_colors = {[.3 .3 1],[1 .3 .3]};
        
        plot_celltypes = any([wide,narrow],1);
        
        figure('Position',[94,122,230,820],'Name',['Unit Depths Site ' num2str(site)]);
        plotSpread(templateDepths - (3840-insertion_depth),'categoryIdx', ...
            celltypes,'categoryColors',use_colors(plot_celltypes));
        set(gca,'XTick',[]);
        set(gca,'YDir','reverse');
        ylabel('Depth (\mum)');
        title(['units by depth, site ' num2str(site)])
        ylim([-100,insertion_depth])
%         plot([-100 100],[max_depth_to_analyze max_depth_to_analyze],'linewidth',2,'linestyle',':','color',[.5 .5 .5])
        
        legend(celltype_labels(plot_celltypes),'location','southoutside');
        set(gca,'Color','k')
        drawnow;
        
    end
    
end


%% Create a cell array of spike times for each unit on each probe, for both good units and good/MUA, then save

% initialize cell array
curr_site_ephys_data = cell(size(good_templates_idx,2));

% populate with spike times for each unit (add LFP, eventually)
for unitID = 1:size(good_templates_idx,1)
    curr_site_ephys_data{unitID,1} = spike_times(spike_clusters==unitID);
    curr_site_ephys_data{unitID,2} = spike_times_timeline(spike_clusters==unitID);
end
ephys_data{site,1} = curr_site_ephys_data;

% add in unitIDs, widths (binary), and depth (relative to brain surface) for the site
ephys_data{site,2} = spike_times;
ephys_data{site,3} = spike_times_timeline;
ephys_data{site,4} = spikeDepths;
ephys_data{site,5} = spike_clusters;
ephys_data{site,6} = good_templates_idx;
ephys_data{site,7} = wide;
ephys_data{site,8} = templateDepths - (3840-insertion_depth);  



end


%% Save and finish

if include_MUA; mua_text = 'mua_';
else; mua_text = '';
end

if max_depth_to_analyze; depth_text = [num2str(max_depth_to_analyze) 'um'];
else; depth_text = '';
end

if save_data
    save([save_folder 'ephys_' mua_text animal '_' day '_' experiments{experiment} '_' depth_text],'ephys_data')
end

if verbose; disp('Finished loading experiment.'); end