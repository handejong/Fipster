%% Automated workflow for making peri-event plots of FIP data


%% You might want to clear workspace
close all
clear
clc


%% Variables
fip_filename='A848_0_000.mat';
data_filename='A848_0';
do_you_want_to_see_FIP_signal=false; % or true
stamps_name='licks';
fiber=1;


%% Getting time stamps, uncomment whatever is applicable

%% Stamps from MED_PC file

mouse=import_medPC(data_filename);
stamps=mouse(end).G; % If stamps are in variable G
stamps=stamps(stamps~=0); % Remove stamps at 0

%.... Maybe some more processing? For instance, only use stamps if there
% was no stamps for 5s before (could be start of bout).
new_stamps=stamps(1); %there was nothing in the 5s before the first stamp
for i=2:length(stamps)
    if stamps(i)-stamps(i-1)>5
        new_stamps=[new_stamps, stamps(i)];
    end
end
stamps=new_stamps;

%% Stamps from somewhere else
% stamps=... (put variable name)


%% Import FIP data
if do_you_want_to_see_FIP_signal
    signal=FIP_signal('Filename',fip_filename);
else
    signal=FIP_signal('Filename',fip_filename,'no figure');
end

%% Smooth 405nm signal
signal.settings.smooth_405=10; %10 datapoints (1 sec at 10Hz sampling)


%% Let the user choose the correct channel for the time offset
% If you constantly use the same channel, just set the variable
% offset_channel. Be carefull however, FIP_signal imports all logAI signals
% it can find, and this number changes sometimes. The reason is that the NI
% Daq sometimes records signals on inputs that are not actually connected.
% (ghosting). If more or less signals are imported into FIP_signal, the
% number of the channel you are looking for might change.
temp_fig=figure('Position',[1,1,800,200]);
for i=2:size(signal.logAI,2)
    plot(signal.logAI(:,1),signal.logAI(:,i),'DisplayName',num2str(i));
    hold on
end
legend
offset_channel=input('What number is the offset channel? ');
delete(temp_fig);

%Or set the variable offset_channel manually:
% offset_channel=4;


%% Derive stamps and import them
scr.YData=signal.logAI(:,offset_channel)';
signal.import_timestamps(signal.derive_stamps(scr),stamps_name);


%% Use stamps to offset signal
stamp_number=1; % Which stamp should be used?
is_time=15; % This stamps is how many seconds?
time_offset=signal.timestamps{end}(stamp_number)-is_time; %
disp(['Time offset is ' num2str(time_offset) ' sec.'])
signal.settings.time_offset=time_offset;


%% Import external stamps
signal.import_timestamps(stamps,stamps_name);
PE_plot=signal.peri_event_plot(signal.data{fiber},signal.timestamps{end},15); %15 is the window size in seconds


%% Set the baseline in the sweepset object and substract
% Note that this object is in ms, while FIP_signal was in seconds.
PE_plot.settings.baseline_info.start=-5000; % -5000ms
PE_plot.settings.baseline_info.end=0; % to 0ms
PE_plot.settings.baseline_info.substracted=true;

% Maybe plot Z-scores instead
PE_plot.settings.Z_scores=true;
PE_plot.refocus();

% Smooth data
PE_plot.smooth_trace(1);  %In ms, but be carefull not to go below the sampling interval (100ms at 10hz sampling) becaue the smooth factor (or smooth span) will be interpreted as a percentage.


%% Plot Heatplot (or whatever that is called)
PE_plot.heat_plot;
