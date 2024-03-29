%% Automated workflow for making peri-event plots of FIP data
% FIPSTER is a GUI based MATLAB interface for fiber-photometry signals
% obtained using multi-fiber camera-based systems. However, it is also
% possible to completely automated FIP analysis using FIPSTER. This file is
% meant to give some pointers on how to do this.

%% You might want to clear workspace
close all
clear
clc


%% Input data
% As example input we have provided an 3-fiber 1-channel example signal
% recording. There are to logAI traces in the example data as well. Here
% we will use them to obtain task-event and make peri-event plots.
fip_filename = 'raw_data/example_1.mat';


%% Settings

% Do you want to see the GUI, or run everyhing in the background? We
% recomend you set this to "true" the first time. This will allow you to
% use the GUI as well as the automated aproach.
do_you_want_to_see_FIP_signal=true;


%% Import FIP data
% This will make the FIP_signal object. You can navigate the settings of
% this object using the GUI, using dot-notation or using the MATLAB object
% inpsector by double-clicking ont he object name in the sidebar on the
% right.
%
% NOTE: because there is no 405nm in the example signal, FIPSTER will ask
% if you want to detrend the signal. It does so by fitting and subtracting
% the signal in the file "mean_bleaching.mat".

if do_you_want_to_see_FIP_signal
    signal=FIP_signal('Filename',fip_filename);
else
    signal=FIP_signal('Filename',fip_filename,'no figure');
end


%% Smooth 405nm signal
% The example data does not have a 405nm (ref) signal, but if it does, you
% can smooth it here if you want that. For instance, smoothing 10
% datapoints (1sec at a framerate of 10Hz) looks like this:

signal.settings.smooth_405=10;


%% Taking a look at the Log_AI signals and TTL pulses
% The logAI traces record everything that is connected to the DAQ system
% during a fiber photometry recording. You can use this to record animal
% behavior for instance. The example data has two log_AI traces. Channel
% one looks at mouse licking behavior and channel 2 looks at the trial
% structure. In channel 2, trial onset (of a Pavlovian task) is signalled
% by the offset of a 1-sec pulse. Reward delivery is signalled by a 100ms
% pulse and CS- onset is signalled by a 0.5-sec pulse.

% You can convert the log_AI traces into timestamps simply by
% right-clicking them in the GUI. Or you can automate it using the buil-in
% method "derive_stamps".
log_AI_1 = signal.handles.logAI_plots{1};
licks = signal.derive_stamps(log_AI_1);

% Note how 'licks' is just an array with a bunch of timestamps. We can are
% now going to import these timestamps into the FIP_signal object. You can
% do this with any similar (externally generated) timestamps even if they
% were not recorded using the DAQ (as long as the timelines are alligned).
signal.import_timestamps(licks, 'licks')

% To look at the timestamps in the GUI, press the "Time stamps" button on
% the left of the logAI (bottom) plot.

% For the trial onset, deriving the timestamps is a little bit more
% complicated because there are stamps of different length in the signal
% and we want to look at the stamps offset, not the onset.
log_AI_2 = signal.handles.logAI_plots{2};
min_pulse_length = 0.75; %s
max_pulse_length = 1.5; %s
rising = false;
trial_onset = signal.derive_stamps(log_AI_2, min_pulse_length,...
                    max_pulse_length, rising);
signal.import_timestamps(trial_onset, 'trial_onset')

% In the same way we can import reward delivery delivery
min_pulse_length = 0.001; %s
max_pulse_length = 0.2; %s
rising = false;
reward = signal.derive_stamps(log_AI_2, min_pulse_length,...
                    max_pulse_length, rising);
signal.import_timestamps(reward, 'reward')

% And CS- offset
min_pulse_length = 0.25; %s
max_pulse_length = 0.75; %s
rising = false;
CS_minus = signal.derive_stamps(log_AI_2, min_pulse_length,...
                    max_pulse_length, rising);
signal.import_timestamps(CS_minus, 'CS-')


%% Timeline offset
% Sometimes you might want to offset the FIP signal to some
% externally-generated TTL pulse, for instance if you use sync pulses to
% aligh different datasets. You can do this using the GUI (right-click and
% use the "use for time alignment. You can also do it by using the setting
% 'time_offset'.

% For instance, imagine that we want the first trial onset to be at T=0.

%first_trial = signal.timestamps{2}(1); <- UNCOMMENT THIS
%signal.settings.time_offset=first_trial; <- UNCOMMENT THIS

% NOTE in this example the timestamps did not move, so you'd have to
% re-derive them (simply by running the above code cell again).


%% Making peri-event plots
% This can be easily done using the GUI (righ click on a timestamp
% sequence) or using code as follows:

% Make a peri-event plot for the trial onset
data = signal.data;
stamps = signal.timestamps{2};
window = 15; %s
PE_plot=signal.peri_event_plot(data, stamps ,15); 


%% Work with the peri-event (sweepset) object
% Code for the peri-event object is in the file "sweepset.m". Note that it
% was orignally made for ex-vivo electrophysiology data and that is why the
% timeline is in ms (not s).

% All the opperations below can be performed by righ clicking on the
% peri-event plot, either on the background or on individual sweeps. Use
% SPACE to switch between different channels (there are 3 in the example
% data).

% Set a baseline from -10 to -2 sec
PE_plot.settings.baseline_info.start=-10000; %ms
PE_plot.settings.baseline_info.end=-2000; % to 0ms
PE_plot.settings.baseline_info.substracted=true;

% Maybe plot Z-scores instead
PE_plot.settings.Z_scores=true;
PE_plot.refocus();

% Smooth data
% NOTE we actually NEVER smooth data. It's unclear why anybody would ever
% want to do this other than to make something "look good".
smooth_span = 1000; %ms NOTE don't go bellow the sampling interval
%PE_plot.smooth_trace(smooth_span);  %<-- UNCOMMENT THIS


%% Finally, you can plot a a Heatplot
for channel = 1:3
    PE_plot.current_channel=channel;
    PE_plot.heat_plot;
end


%% Or output the data for use somewhere else
% Try 'whole_trace' instead of average. The output will be assigned in the
% workspace. The first column is always a timeline.
for channel = 1:3
    PE_plot.current_channel = channel;
    name = ['channel_', num2str(channel)];
    PE_plot.output_data('average', name)
end
