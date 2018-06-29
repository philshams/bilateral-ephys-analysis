
function [firing_rates, trial_tertiles] = get_avg_frs(spikeTimes, eventTimes, window, num_neurons);   

%
% Fast computation of avg firing rate and spike counts in a window relative to some
% events. 
%


spikeTimes = spikeTimes(:);
eventTimes = sort(eventTimes(:));

% first we'll subselect spikes that are between the end and beginning of
% all the ranges - in some cases this will be useless but often instead
% reduces spike count by a lot.
spikeTimes = spikeTimes(spikeTimes>min(eventTimes+window(1)) & spikeTimes<max(eventTimes+window(2)));
binsize = window(2) - window(1);
[binnedArray, bins] = timestampsToBinned(spikeTimes, eventTimes, binsize, window);
firing_rates = (binnedArray./binsize) / num_neurons;

% get array of 1, 2, or 3 corresponding to firing rate tertiles
trial_tertiles = zeros(length(firing_rates),1);
tertile_thresholds = prctile(firing_rates,[33 66]);
trial_tertiles(firing_rates <= tertile_thresholds(1)) = 1;
trial_tertiles(firing_rates > tertile_thresholds(1) & firing_rates <= tertile_thresholds(2)) = 2;
trial_tertiles(firing_rates > tertile_thresholds(2)) = 3;









