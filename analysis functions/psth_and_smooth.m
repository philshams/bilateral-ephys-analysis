

function [psth_smooth_all, psth_smooth, psth_std_trials, bins, rasterX, rasterY, spikeCounts, binnedArray] = ...
                                                        psth_and_smooth(spikeTimes, eventTimes, window, psthBinSize, smooth_window, num_neurons)
% function [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikeTimes, eventTimes, window, psthBinSize)
%
% Fast computation of psth and spike counts in a window relative to some
% events. Also returns rasters you can plot. 
%
% Notes on inputs:
% - eventTimes is nEvents x 1
% - window is length 2, e.g. [-0.1 0.3] for a window from -0.1 to +0.3 sec
% relative to events. 
%
% Notes on outputs:
% - psth can be plotted with plot(bins, psth);
% - rasters can be plotted with plot(rasterX, rasterY);
% - spikeCounts is nEvents x 1, where each entry is the number of spikes
% that occurred within the window around that event. 

spikeTimes = spikeTimes(:);
eventTimes = sort(eventTimes(:));

% first we'll subselect spikes that are between the end and beginning of
% all the ranges - in some cases this will be useless but often instead
% reduces spike count by a lot.
spikeTimes = spikeTimes(spikeTimes>min(eventTimes+window(1)) & spikeTimes<max(eventTimes+window(2)));
[binnedArray, bins] = timestampsToBinned(spikeTimes, eventTimes, psthBinSize, window);
firing_rates = (binnedArray./psthBinSize) / num_neurons;

% make psth
psth = mean(firing_rates, 1); % normalize to Hz

% smooth psth by trial to get std
if smooth_window
    gw = gausswin(smooth_window)';
    smWin = gw./sum(gw);
    psth_smooth_all = conv2(firing_rates, smWin, 'same');
else
    psth_smooth_all = firing_rates;
end
psth_smooth = mean(psth_smooth_all,1);


% get std across trials
if size(psth_smooth_all,1) == 1
    psth_std_trials = 0;
else
    psth_std_trials = std(psth_smooth_all, 1);
end

% addtl saved outputs...
[tr,b] = find(binnedArray);
[rasterX,yy] = rasterize(bins(b));
rasterY = yy+reshape(repmat(tr',3,1),1,length(tr)*3); % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything
spikeCounts = sum(binnedArray,2);

% customize psth etc for plotting purposes
psth_std_trials = psth_std_trials(:,round(smooth_window/2)+1:end-round(smooth_window/2))';
psth_smooth = psth_smooth(:,round(smooth_window/2)+1:end-round(smooth_window/2))';
psth_smooth_all = psth_smooth_all(:,round(smooth_window/2)+1:end-round(smooth_window/2));
bins = bins(round(smooth_window/2)+1:end-round(smooth_window/2))';


