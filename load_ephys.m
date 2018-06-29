% -----------------------------------------------------------------
%                 Load and Examine Ephys Data
% -----------------------------------------------------------------


%% LOAD DATA

save_folder = 'P:\ephys data\';
save_data = true;

animal = 'DS_PS001'; % animal name in string
day = '2018-05-15'; % date in string
experiment = 3; % experiment number (sparse noise, full-screen flicker, choice world, gray screen)
experiments = {'noise','flicker','choice','spontaneous'};
verbose = true; % prints progress of load
include_MUA = false;
insertion_depth = 2000;
max_depth_to_analyze = 1000;

sites = 1:2;
site_sides = {'right','left'};






%% PROCESS DATA

% load experiment details
AP_load_experiment;

% only need to run this once per experiment / include_MUA / max_depth_to_analyze setting -- can load using the next cell thereafter
% load_save_ephys;


%% LOAD PROCESSED DATA

% include MUA? limit depth?
if include_MUA; mua_text = 'mua_';
else; mua_text = ''; end

if max_depth_to_analyze; depth_text = ['_' num2str(max_depth_to_analyze) 'um'];
else; depth_text = ''; end

% load data
ephys_data = load([save_folder 'ephys_' mua_text animal '_' day '_' experiments{experiment} depth_text]);
disp([ 'ephys_' mua_text animal '_' day '_' experiments{experiment} depth_text ' loaded'])




%% ANALYZE DATA
close all

epoch_to_split_by = [0 .5]; % time in s relative to stimulus to examine contralateral spontaneous activity
smooth_window = 100; % size of gaussian window for PSTH


% ---------  plot trials (stim type and each hemisphere's z-scored-by-stimulus activity third) -----------
show_behaviour;


% ----------- plot population PSTH -----------
stims_by_side = {[1,2,3,4,5], [6,7,8,9,10]}; % remove the 1 and 6 to get rid of so-called 0-contrast trials
contrasts = {'0%','12%','25%','50%','100%'}; % remove the '0%' to get rid of so-called 0-contrast trials
stim_sides = {'left','right'};

plot_correlations = false; % also plot correlations over time -- I didn't find this particularly enlightening

separate_contralateral_spontaneous = false;
stimulus_response;   
separate_contralateral_spontaneous = true;
stimulus_response;   


% ----------- plot single unit PSTH ----------- saves figures to the 'save_folder'
single_unit_stimulus_response;


% ----------- plot contrast response function -----------
epoch_to_split_by = [0 .5]; 
contrast_response_function % for all trials
contrasts_response_by_contra % separate by thirds

epoch_to_split_by = [0 .2]; 
contrast_response_function % for all trials
contrasts_response_by_contra % separate by thirds


%%  ----------- look at correlations and variance explained -----------

% parameters
shift_back_time = 0; % allows you to look at correlations offset in time
stims_by_side =  {[7,8,9,10,2,3,4,5], [2,3,4,5,7,8,9,10]};
title_suffix = '';

variance_by_baseline = true; % use coefficient of contra activity using the baseline period
analyze_single_units = false;

baseline_correlation % get the coefficients to use to predict activity, and get correlations during baseline period


% get avg baseline variance explained by stimulus to subtract from what you
% get during the stimulus; and same for contra activity with shuffled data,
% in the case of NOT using the coefficient derived from the baseline. This takes a
% while if analyze_single_units = true;
stims_by_side =  {[7,8,9,10,2,3,4,5], [2,3,4,5,7,8,9,10]};
variance_explained_get_baseline;

% get the results
variance_explained_over_time;


% do so for weak stimuli
stims_by_side = {[7,2], [2,7]};
title_suffix = '- weak stimuli';
baseline_correlation;
variance_explained_get_baseline;
variance_explained_over_time;

% do so for strong stimuli
stims_by_side = {[10,5], [5,10]};
title_suffix = '- strong stimuli';
baseline_correlation;
variance_explained_get_baseline;
variance_explained_over_time;



%%  ------------ do single unit analyses (some steps here take a while)  ------------

% parameters
stims_by_side =  {[7,8,9,10,2,3,4,5], [2,3,4,5,7,8,9,10]};
title_suffix = '';
variance_by_baseline = true; % use coefficient of contra activity using the baseline period
analyze_single_units = true;

% get the coefficient etc. for each unit
baseline_correlation; 
stims_by_side = {[2,3,4,5], [7,8,9,10]}; 
baseline_correlation_stimulus; % get the coefficient etc. for each unit during stimuli

stims_by_side =  {[7,8,9,10,2,3,4,5], [2,3,4,5,7,8,9,10]};
variance_explained_get_baseline; % get the baseline variance explained -- takes a while


% scatter plot for each neuron and population of variance explained
epoch_to_split_by = [0 .5]; 
variance_explained_fr; 

% get data on units' firing rates and responsiveness for plots below
stims_by_side = {[2,3,4,5], [7,8,9,10]}; 
response_by_depth % adds info to baseline_correlation for each neuron (firing rate, responsiveness) to be used in correlation_by_depth



% correlation by depth plot
during_stim = false; % use baseline_correlation_stimulus or baseline_correlation;
show_responsiveness = false; %circle vs square
show_baseline_fr = true; % size of point

correlation_by_depth

during_stim = true; % also show during stimulus
correlation_by_depth


% look at the effect of stimulus on correlations
stims_by_side = {[2,3], [7,8]};
baseline_correlation_stimulus_weak; 
stims_by_side = {[4,5], [9,10]};
baseline_correlation_stimulus_strong;  
% NOTE: for the last experiment, upper depths were more modulated during strong stimuli, and lower depths more modulated during weak stimuli. 
    % It would be interesting to see if this pattern persists.

use_stim_type = 'all';
correlation_by_depth_stim_effect
use_stim_type= 'weak';
correlation_by_depth_stim_effect
use_stim_type = 'strong';
correlation_by_depth_stim_effect



%%  some extra measures of correlation!

stims_by_side = {[2,3,4,5], [7,8,9,10]};
epoch_to_correlate = [0 .5];
show_correlation_plot = true;

depth_correlation
contralateral_depth_correlation

unit_correlation
contralateral_unit_correlation

pairwise_correlations_over_time % this one takes a while























