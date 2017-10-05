classdef FIP_signal <handle
    %FIP_signal Imports fiber photometry signal into FIPster
    %   FIP_signal takes either a .mat file or one or two variables as
    %   input. It will attempt to contruct valid fiber photometry signals
    %   from the input for multiple fibers and up to three channels. If one
    %   of the data channels is a signal collected after emission with
    %   405nm light, it will use that signal to correct the other channels
    %   for movement artifacts.
    %
    %
    %   INPUT ARGUMENTS:
    %       - 'Filename',filename       filename should refer to a .mat
    %                                   file which will be imported.
    %       - 'User input'              Will open a file browser window.
    %       - 'Signal',sig              Sig should be and ix1 array were
    %                                   i represents the # of datapoints.
    %                                   Multiple signals can be imported in
    %                                   this way.                               
    %       - 'Refference',ref          Ref should be formated as signal.
    %                                   Note that refference traces (405nm 
    %                                   signal) should be presented
    %                                   in the same order as signal traces.
    %       - 'Time',time_info          Time info can be either a complete
    %                                   timeline formatted ix1, one value,
    %                                   which will be asumed to be the
    %                                   sampling frequency or 2 values wich
    %                                   will be presumed to be the start
    %                                   (first timepoint) and end (last
    %                                   timepoint) of the associated
    %                                   signal. The user can either present
    %                                   one time signal which will be
    %                                   applied to all signals or one time
    %                                   argument for each signal, but no
    %                                   other number of time arguments.
    %       - 'Signal set',signal       Signal is an i,2,n matrix were i is
    %                                   datapoints and n is the number of
    %                                   fibers. :,1,: should contain signal
    %                                   and :,2,: should contain a
    %                                   timeline. If signal is formatted
    %                                   i,3,n, :,1: will be assumed to be a
    %                                   raw signal, :,2,: to be a
    %                                   refference signal and :,3: to be a
    %                                   timeline.
    %       - 'no figure'               The object is still created, but no
    %                                   figure. This is efficient for other
    %                                   programs (e.g. Fipster).
    %
    %   Note: 'Filename','User_input' and 'Signal_set' are not compattible
    %   with 'Signal, 'Refference' and 'Time'.
    %
    %
    %   METHODS
    %       - 'go_to_time'              Will scroll to the nearest available
    %                                   timepoints or interpolate when this 
    %                                   is turned on.
    %
    %
    % FIP_signal is part of FIPster. FIPster is made by Johannes de Jong,
    % j.w.jong@berkeley.edu
     
    properties
        info                    % Info about the current data and raw data
        raw_data                % Contains data as imported on contruction
        timestamps              % Lists of timestamps (imported or generated)
        timestamps_names        % Names of the timestamp traces
        np_signals              % Number of signals
        handles                 % Of all figure ellements
        window                  % 'all time' or time window in sec.
    end
    
    properties (Dependent)
        data                    % Ouput data ix2
        sig_CD                  % Calcium-dependend signals
        sig_405                 % Non-CD signals
        logAI                   % Logdata
    end
    
    properties (Access = private)
        mouse_click             % Boolean if the mouse is pressed
        mouse_start             % Mouse start coordinates
        mouse_action            % Struct with mouse action info
        r_mouse_scr             % Source of a right mouse click (context)
        figure_open             % Boolean stores if there is a figure
        raw_signal              % Boolean if raw signal or data is shown
        AI_plots                % Boolean if shown analog input or timestamps
        fibers_shown            % Boolean array, which fibers displayed
        log_AI_used             % Boolean array, which logAI not deleted
        key_down                % If a key is down on the main figure
        key_pressed             % Pressed key
    end
    
    properties (SetObservable)
        settings                % Contains all settings
        c_time                  % Current time
        c_signal                % Current signal
    end
    
    methods
        
        function this_FIP_signal=FIP_signal(varargin)
            % Constructor
            
            % Default variables (can be changed by input arguments)
            present_figure=true;
            this_FIP_signal.mouse_click=false;
            this_FIP_signal.raw_signal=true;
            this_FIP_signal.AI_plots=true;
            this_FIP_signal.c_time=0;
            
            % Populate settings (can be changed by input arguments)
            settings.time_offset=0;
            settings.smooth_405=1;
            settings.smooth_CD=1;
            settings.smooth_data=1;
            settings.fit_405='polyfit_2';
            this_FIP_signal.settings=settings;
            
            % Dealing with input arguments
            skipp_next=false;
            data_loaded=false;
            for i=1:nargin
                if ~skipp_next
                    switch varargin{i}
                        case 'Filename'
                            if data_loaded
                                error('Please choose only one method to import data, either User input, Filename or individual variables')
                            end
                            path=which(varargin{i+1});
                            path=path(1:end-length(varargin{i+1}));      
                            this_FIP_signal.load_file(varargin{i+1},path);
                            skipp_next=true;
                            data_loaded=true;
                        case 'User input'
                            if data_loaded
                                error('Please choose only one method to import data, either User input, Filename or individual variables')
                            end
                            [filename, path]=uigetfile({'*.mat'},'select file','MultiSelect','off');
                            this_FIP_signal.load_file(filename,path);
                            data_loaded=true;
                        case 'no figure'
                            present_figure=false;
                        otherwise
                            warning(['input ' varargin{i} ' not recognized.'])
                    end
                else
                    skipp_next=false;
                end
            end
            disp(['Loaded ' num2str(this_FIP_signal.np_signals) ' signals.'])
            
            % Error handling <- but not the errors in load_data.... TODO
            if ~data_loaded
                error(['Please choose a valid method to import data, either User input, Filename or individual variables.' ...
                '(FIP_signal is case sensitive.)'])
            end
                    
            % Present figure
            if present_figure
                this_FIP_signal.present_figure();
            end    
        end
        
        function load_file(this_FIP_signal, filename, path)
            % Will load a .mat file into the FIP_signal object
            this_FIP_signal.np_signals=0;
            
            if ~strcmp(filename(end-3:end),'.mat')
                error 'please submit a .mat file'
            end
            
            load([path filename])
            
            variables=who;
            
            for i=1:length(variables)
                
                switch variables{i}
                    case 'sig'
                        for j=1:length(sig(1,:))
                            this_FIP_signal.raw_data.sig{this_FIP_signal.np_signals+1}=sig(1:end-5,j);
                            this_FIP_signal.raw_data.ref{this_FIP_signal.np_signals+1}=ref(1:end-5,j);
                            this_FIP_signal.np_signals=this_FIP_signal.np_signals+1;
                        end
                    case 'framerate'
                        this_FIP_signal.info.framerate=framerate/2;
                    case 'labels'
                        for j=1:length(labels)
                            this_FIP_signal.info.names{j}=labels{j};
                        end
                    case {'this_FIP_signal' 'filename' 'path' 'ref'}
                        % Good, but do nothing
                        all_good=true;
                    otherwise
                        warning(['Variable ' variables{i} ' was not imported'])
                end     
            end
            
            % Do we find any associated files on the same path
            files=dir(path);
            for i=1:length(files) 
                if length(files(i).name)>length(filename) && strcmp(files(i).name(1:length(filename)-4),filename(1:end-4))
                    if length(files(i).name)>8 && strcmp(files(i).name(end-8:end),'logAI.csv')
                        % found a logAI file, check if not to big!
                        logAI=csvread([path files(i).name]);
                        disp(['Loaded ' files(i).name])
                        logAI=logAI(:,max(logAI)>1); %only import channels with signal >1mV
                        this_FIP_signal.raw_data.logAI=logAI;
                        dims=size(logAI);
                        disp(['Loaded ' num2str(dims(2)) ' analog inputs.']);
                        this_FIP_signal.log_AI_used=true(dims(2),1);
                        this_FIP_signal.AI_plots=true;
                    end
                    if length(files(i).name)>14 && strcmp(files(i).name(end-14:end),'calibration.jpg')
                        % found a callibration picture
                        this_FIP_signal.info.callibration=imread([path files(i).name]);
                        disp(['Loaded ' files(i).name])
                    end
                end
            end
            
            % Populate object info (Note: of the raw data!)
            for i=1:this_FIP_signal.np_signals
                this_FIP_signal.info.signal_units{i,1}='F';    %for signals    
                this_FIP_signal.info.signal_units{i,2}='deltaF/F'; %for data
                this_FIP_signal.info.max_time(i)=length(this_FIP_signal.raw_data.sig{i}(1:end-1))/this_FIP_signal.info.framerate;
            end
            
            % Populate settings
            m_settings=this_FIP_signal.settings;
            m_settings.signal_units=this_FIP_signal.info.signal_units;
            this_FIP_signal.settings=m_settings;
        end
        
        function present_figure(this_FIP_signal)
            % draws a figure to display raw and precessed data. Also writes
            % all the context menus.
            
            % Present main figure
            this_FIP_signal.handles.main_figure=figure('name','FIP signals',...
                'WindowKeyPressFcn',@this_FIP_signal.key_press,...
                'WindowKeyReleaseFcn',@this_FIP_signal.key_release,...
                'Tag','main_figure',...
                'Position',[1 1 1317 689],...
                'CloseRequestFcn',@this_FIP_signal.close_req,...
                'Menubar','none');
            
            % Populate context (drop-down) menus
            this_FIP_signal.make_context();
            
            % Populate the fibers_shown boolean
            this_FIP_signal.fibers_shown=false(this_FIP_signal.np_signals,1);
            
            % Plotting of the signals
            names=this_FIP_signal.info.names;
            this_FIP_signal.handles.FIP_plot=subplot(...
                'Position',[0.13 0.3 0.835 0.65],...
                'ButtonDownFcn',@this_FIP_signal.mouse_down,...
                'Tag','FIP_plot',...
                'UIContextMenu',this_FIP_signal.handles.drop_down.FIP_plot);
            hold on
            for i=1:this_FIP_signal.np_signals
                this_FIP_signal.handles.sig_CD{i}=plot(this_FIP_signal.sig_CD{i}(:,2),this_FIP_signal.sig_CD{i}(:,1),...
                    'Visible','off',...
                    'DisplayName',[names{i} ' raw signal'],...
                    'Color',[0 0.45 0.75],...
                    'ButtonDownFcn',@this_FIP_signal.mouse_down,...
                    'Tag',['sig_CD.' num2str(i)]);
                hold on
                this_FIP_signal.handles.sig_405{i}=plot(this_FIP_signal.sig_405{i}(:,2),this_FIP_signal.sig_405{i}(:,1),...
                    'Visible','off',...
                    'DisplayName',[names{i} ' 405 signal'],...
                    'Color',[0.85 0.32 0.1],...
                    'ButtonDownFcn',@this_FIP_signal.mouse_down,...
                    'UIContextMenu',this_FIP_signal.handles.drop_down.sig_405,...
                    'Tag',['sig_405.' num2str(i)]);
                this_FIP_signal.handles.data{i}=plot(this_FIP_signal.data{i}(:,2),this_FIP_signal.data{i}(:,1),...
                    'Visible','off',...
                    'DisplayName',[names{i} ' corrected'],...
                    'ButtonDownFcn',@this_FIP_signal.mouse_down,...
                    'UIContextMenu',this_FIP_signal.handles.drop_down.data,...
                    'Tag',['data.' num2str(i)]);
            end
            this_FIP_signal.c_signal=1;
            this_FIP_signal.handles.sig_CD{this_FIP_signal.c_signal}.Visible='on';
            this_FIP_signal.handles.sig_405{this_FIP_signal.c_signal}.Visible='on';
            this_FIP_signal.handles.ylabel=ylabel(this_FIP_signal.settings.signal_units{1,1});
            xlabel('Time (s)');
            
            % Plotting the line indicating the current time
            temp=ylim;
            m_c_time=this_FIP_signal.c_time;
            
            this_FIP_signal.handles.c_time_line=line([m_c_time m_c_time],[temp(1) temp(2)],...
                'LineStyle','--',...
                'Color',[0.2 0.2 0.2],...
                'Visible','off');
           
            % logAI figure
            this_FIP_signal.handles.logAI_plot=subplot(...
                'Position',[0.13 0.04 0.835 0.18],...
                'Tag','logAI_plot',...
                'ButtonDownFcn',@this_FIP_signal.mouse_down,...
                'UIContextMenu',this_FIP_signal.handles.drop_down.FIP_plot);
            hold on
            
            % Plotting logAI signals
            logAI_plot=this_FIP_signal.raw_data.logAI;
            dims=size(logAI_plot);
            for i=2:dims(2) % first column is the timeline
                if this_FIP_signal.log_AI_used(i)
                    this_FIP_signal.handles.logAI_plots{i}=plot(logAI_plot(:,1),logAI_plot(:,i),...
                        'DisplayName',['AI ' num2str(i)],...
                        'UIContextMenu', this_FIP_signal.handles.drop_down.logAI,...
                        'ButtonDownFcn',@this_FIP_signal.mouse_down,...
                        'Tag',['AI ' num2str(i)]);
                else
                    this_FIP_signal.handles.logAI_plots{i}='deleted logAI plot, see raw data';
                end
            end
            ylabel('(mV)')
            
            % Plotting any timestamp data (if available)
            if isfield(this_FIP_signal.handles,'time_stamp_plots')
                for i=1:length(this_FIP_signal.timestamps)
                    this_FIP_signal.handles.time_stamp_plots{i}=...
                        plot(this_FIP_signal.timestamps{i},...
                        ones(length(this_FIP_signal.timestamps{i}),1)*i,'x');
                end
            end
            
            
            % Text to indicate which fiber or data the user is looking at
            this_FIP_signal.handles.text_1=uicontrol(gcf,...
                'Units','normalized',...
                'Position',[0.88 0.9 0.05 0.05],...
                'String',this_FIP_signal.info.names{1},...
                'Style','text',...
                'FontSize',20,...
                'BackgroundColor',[1 1 1]);
            
            % Let's add some buttons
            uicontrol(this_FIP_signal.handles.main_figure,...
                'Units','normalized',...
                'Position',[0.01 0.9400 0.0583 0.0301],...
                'Style','text',...
                'String','Display:',...
                'HorizontalAlignment','left');
            
            % Fiber display selection (checkboxes)
            for i=1:this_FIP_signal.np_signals % for every fiber
                j=i-1;
                uicontrol(this_FIP_signal.handles.main_figure,...
                'Units','normalized',...
                'Position',[0.01 0.91-(j*0.03) 0.0583 0.0301],...
                'Style','checkbox',...
                'String',this_FIP_signal.info.names{i},...
                'Tag',['signal_' num2str(i)],...
                'Callback',@this_FIP_signal.UI_used,...
                'HorizontalAlignment','left');
            end
            
            % Timestamps or logAI button
            this_FIP_signal.handles.time_stamp_button=uicontrol(...
                this_FIP_signal.handles.main_figure,...
                'Units','normalized',...
                'Position',[0.01 0.1899 0.0583 0.0301],...
                'Style','togglebutton',...
                'String','Time stamps',...
                'Tag','stamps_switch',...
                'Callback',@this_FIP_signal.UI_used,...
                'HorizontalAlignment','left');
            
            % Slider for time control
            max_time=max(this_FIP_signal.info.max_time);
            sliderstep=[1/max_time 10/max_time];
            this_FIP_signal.handles.time_bar=...
                uicontrol(this_FIP_signal.handles.main_figure,...
                'Units','normalized',...
                'Position',[0.13 0.22 0.835 0.029],...
                'Style','slider',...
                'Tag','time_slider',...
                'Callback',@this_FIP_signal.UI_used,...
                'Max',max_time,...
                'SliderStep',sliderstep);
            
            % norm. or raw data button
            this_FIP_signal.handles.norm_data_button=...
                uicontrol(this_FIP_signal.handles.main_figure,...
                'Units','normalized',...
                'Position',[0.01 0.91-(i*0.03)-0.02 0.0583 0.0301],...
                'Style','togglebutton',...
                'String','norm. data',...
                'Tag','data_switch',...
                'Callback',@this_FIP_signal.UI_used,...
                'HorizontalAlignment','left');

            % Add listeners
            this_FIP_signal.handles.listeners{1}=addlistener(this_FIP_signal,...
                'settings','PostSet',@this_FIP_signal.update_plots);    
            
            % Wrapping up
            this_FIP_signal.figure_open=true;
            this_FIP_signal.window='all time';
            
            % Immediately deal with dependend variables
            this_FIP_signal.update_plots;
        end
        
        function stamps=derive_stamps(this_FIP_signal, scr)
            % Converts input signal (scr) to time stamps
            
            % What kind of input signal?
            scr_type='AI'; % only possibility for now
            
            % Default variables
            start_pulse = true; % take start or end of pulse
            threshold = 0.1+3*std(scr.YData);
            pulse_per_stamp=1;
            
            % Figure out if the pulses go down or up
            max_signal=max(scr.YData);
            min_signal=min(scr.YData);
            mean_signal=mean(scr.YData);
            
            if abs(mean_signal-min_signal)>abs(mean_signal-max_signal)
                threshold=mean_signal-threshold;
                stamp_logic=scr.YData<threshold;
                disp(['Pulse is defined as signal <' num2str(threshold) ' mV.']);
            else
                threshold=mean_signal+threshold;
                stamp_logic=scr.YData>threshold;
                disp(['Pulse is defined as signal >' num2str(threshold) ' mV.']);
            end

            % Get info from user about conversion and calculate from that
            if strcmp(scr_type,'AI')
                
                % Ask about minium stamp length
                %min_length=inputdlg('
       
                time=this_FIP_signal.logAI(:,1);
                
                
                if start_pulse
                    stamp_logic= stamp_logic & ~[0 stamp_logic(1:end-1)];
                else
                    stamp_logic= ~stamp_logic & [0 stamp_logic(1:end-1)];
                end
                
                stamps=time(stamp_logic,1);    
            end
           

            
        end
        
        function import_timestamps(this_FIP_signal,stamps,name)
            % Will import timestamps and plot them
            
            % Find out how many timestamps arrays there allready are
            nr=length(this_FIP_signal.timestamps)+1;
            
            % Store timestamps and plot them
            this_FIP_signal.timestamps{nr}=stamps;
            this_FIP_signal.timestamps_names{nr}=name;
            subplot(this_FIP_signal.handles.logAI_plot)
            this_FIP_signal.handles.time_stamp_plots{nr}=...
                plot(stamps,ones(length(stamps),1)*nr,'x',...
                'UIContextMenu',this_FIP_signal.handles.drop_down.time_stamps,...
                'Tag',['stamps_' num2str(nr)],...
                'ButtonDownFcn',@this_FIP_signal.mouse_down);
            
            % Done
            disp(['Timestamps ' name ' stored in cell ' num2str(nr) '.'])

            this_FIP_signal.update_plots;
        end
        
        function peri_event_plot(this_FIP_signal, input, stamps, window)
            % Makes a peri-event plot (sweepset) from input signal and
            % stamps (events). Inputs should include a timeline in the 2th
            % column. Time should be in s in stamps,input and
            % window.
            
            % Note: this function assumes equal framerate over the whole
            % input.
            
            % Basic variables
            number_of_events=length(stamps);
            framerate=round(length(input(:,2))/(input(end,2)-input(1,2)));
            timeline=[-window:(1/framerate):window]*1000; %from s to ms
            original_window=window;
            window=window*framerate;
            
            % Check if sufficient data is available around every timestamps
            % (maybe timestamps right at begining or end of session).
            check_finished_start=false;
            check_finished_end=false;
            
            while ~check_finished_start %First early stamps
                if stamps(1)<original_window
                    stamps=stamps(2:end);
                    warning('Removed early time stamp, because no FIP data available.')
                else
                    check_finished_start=true;
                end  
            end
            
            end_time=input(end,2);
            
            while ~check_finished_end
                if (end_time-stamps(end))<original_window
                    stamps=stamps(1:end-1);
                    warning('Removed time stamp at the end, because no FIP data available.')
                else
                    check_finished_end=true;
                end
            end
            
            number_of_events=length(stamps);
            
            results=zeros(number_of_events+1,2*window+1);
            results(1,:)=timeline;
            
            for i=1:number_of_events
                [~, index]=min(abs(stamps(i)-input(:,2)));
                results(i+1,:)=input(index-window:index+window,1);
            end
            test=sweepset('other data',results);
            assignin('base','hoihoi',test)
        end
        
        function sig_CD=get.sig_CD(this_FIP_signal)
            % Get sig_CD based on the settings
            fs=1/(this_FIP_signal.info.framerate);
            for i=1:this_FIP_signal.np_signals
                end_time=(length(this_FIP_signal.raw_data.sig{i})-1)*fs;
                sig_CD{i}(:,1)=smooth(this_FIP_signal.raw_data.sig{i},...
                    this_FIP_signal.settings.smooth_CD);
                sig_CD{i}(:,2)=[0:fs:end_time]';
                
                % Adjust the timeline based on the time_offset.
                sig_CD{i}(:,2)=sig_CD{i}(:,2)-this_FIP_signal.settings.time_offset;
            end
        end
        
        function sig_405=get.sig_405(this_FIP_signal)
            % Gets the sig_405 trace based on the settings
            fs=1/(this_FIP_signal.info.framerate);
            for i=1:this_FIP_signal.np_signals
                end_time=(length(this_FIP_signal.raw_data.ref{i})-1)*fs;
                temp=smooth(this_FIP_signal.raw_data.ref{i},this_FIP_signal.settings.smooth_405);
                switch this_FIP_signal.settings.fit_405
                    case 'unfit'
                        % do not fit 405nm signal
                        sig_405{i}(:,1)=temp;
                        sig_405{i}(:,2)=[0:fs:end_time]';
                    case 'polyfit_2'
                        % default, fit first 2 polynomal coeficients
                        p=polyfit(temp,this_FIP_signal.raw_data.sig{i},1);
                        sig_405{i}(:,1)=temp*p(1)+p(2);
                        sig_405{i}(:,2)=[0:fs:end_time]';
                    case 'polyfit_1'
                        % fit only first polynomal coeficient
                        p=polyfit(temp,this_FIP_signal.raw_data.sig{i},1);           
                        sig_405{i}(:,1)=temp*p(1);
                        sig_405{i}(:,2)=[0:fs:end_time]';
                    case 'fit_mean'
                        % substract the difference of the means
                        sig_mean=mean(this_FIP_signal.raw_data.sig{i});
                        ref_mean=mean(temp);
                        mean_dif=sig_mean-ref_mean;
                        sig_405{i}(:,1)=temp+mean_dif;
                        sig_405{i}(:,2)=[0:fs:end_time]';
                    case 'sliding_window'
                        % polyfit_2, but using a sliding window
                        window_size=400;
                        window_size=window_size-1;
                        sig_405{i}=zeros(length(temp),2);
                        for j=1:window_size:length(temp)-window_size
                            p=polyfit(temp(j:j+window_size),this_FIP_signal.raw_data.sig{i}(j:j+window_size),1);
                            sig_405{i}(j:j+window_size,1)=temp(j:j+window_size)*p(1)+p(2);
                        end
                        j=j+window_size;
                        p=polyfit(temp(j:end),this_FIP_signal.raw_data.sig{i}(j:end),1);
                        sig_405{i}(j:end,1)=temp(j:end)*p(1)+p(2);
                        sig_405{i}(:,2)=[0:fs:end_time]';  
                    otherwise
                        warning('fit method unkown')
                end
                
                % Ajust time based on offset
                sig_405{i}(:,2)=sig_405{i}(:,2)-this_FIP_signal.settings.time_offset;
            end            
        end
        
        function data=get.data(this_FIP_signal)
            % Gets the processed FIP_signal based on the settings
            for i=1:this_FIP_signal.np_signals
                signal_norm(:,1)=this_FIP_signal.sig_CD{i}(:,1)-this_FIP_signal.sig_405{i}(:,1)+this_FIP_signal.sig_CD{i}(1,1);

                % check what unit the data should be, row two is SU for data.
                switch this_FIP_signal.settings.signal_units{i,2}
                    case 'deltaF/F'
                        % Standard F = mean of whole signal
                        % However, later there will be an option to use a
                        % different mean.
                        data{i}(:,1)=(signal_norm(:,1)-mean(signal_norm(:,1)))./mean(signal_norm(:,1));
                    case 'F'
                        % Present F, which is the movement and
                        % bleaching-corrected signal.
                        data{i}(:,1)=signal_norm(:,1);
                    case 'Z-score'
                        % Taking the Z-score directly from the normalized
                        % trace, but I can prove that this is equal to the
                        % Z-score from the deltaF/F trace.
                        data{i}(:,1)=(signal_norm(:,1)-mean(signal_norm(:,1)))/std(signal_norm(:,1));
                end
                
                % Timeline
                data{i}(:,2)=this_FIP_signal.sig_CD{i}(:,2);

                % Apply smooth if requested
                data{i}(:,1)=smooth(data{i}(:,1),this_FIP_signal.settings.smooth_data);
            end       
        end  
        
        function logAI=get.logAI(this_FIP_signal)
            % Get AI log based on settings and selection
            
            % Insert selector boolean before to only load usefull channels
            m_logAI=this_FIP_signal.raw_data.logAI(:,this_FIP_signal.log_AI_used);
            
            % Apply time offset
            m_logAI(:,1)=m_logAI(:,1)-this_FIP_signal.settings.time_offset;
            
            logAI=m_logAI;
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%%%% Callbacks & other functions %%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = private)
        
        function update_plots(this_FIP_signal, varargin)
            % Will update all plots
            
            notify(this_FIP_signal,'state_change');
            
            if ~this_FIP_signal.figure_open % if figure exists (could also just put return)
                % Run getter functions once
                this_FIP_signal.sig_CD;
                this_FIP_signal.sig_405;
                this_FIP_signal.data;
                this_FIP_signal.logAI;
                return
            end
            
            type='standard';
            n_arg=nargin;
            if n_arg>1
                for i=1:(n_arg-1)
                    if ischar(varargin{i}) && strcmp(varargin{i},'type')
                        type=varargin{i+1};
                        break
                    end
                end
            end
            
            switch type % This one goes all the way down
                case 'standard'
                    % Run getter functions once
                    m_sig_CD=this_FIP_signal.sig_CD;
                    m_sig_405=this_FIP_signal.sig_405;
                    m_data=this_FIP_signal.data;
                    logAI_plot=this_FIP_signal.logAI; 

                    % Update plot if the figure exists
                    if sum(this_FIP_signal.fibers_shown)==0
                        use_current=true;
                    else
                        use_current=false;
                    end

                    % Check if all signals have same ylabel of more comples
                    same_y_label.raw=true;
                    same_y_label.data=true;
                    for j=2:this_FIP_signal.np_signals
                        if  ~strcmp(this_FIP_signal.settings.signal_units{j,1}, ...
                                this_FIP_signal.settings.signal_units{j-1,1})
                            same_y_label.raw=false;
                        end
                        if  ~strcmp(this_FIP_signal.settings.signal_units{j,2}, ...
                                this_FIP_signal.settings.signal_units{j-1,2})
                            same_y_label.data=false;
                        end
                    end

                    for i=1:this_FIP_signal.np_signals
                        if (use_current && i==this_FIP_signal.c_signal) || this_FIP_signal.fibers_shown(i)
                            visibility='on';
                        else
                            visibility='off';
                        end

                        if this_FIP_signal.raw_signal
                            set(this_FIP_signal.handles.sig_CD{i},...
                                'Visible',visibility,...
                                'XData',m_sig_CD{i}(:,2),...
                                'YData',m_sig_CD{i}(:,1));
                            set(this_FIP_signal.handles.sig_405{i},...
                                'Visible',visibility,...
                                'XData',m_sig_405{i}(:,2),...
                                'YData',m_sig_405{i}(:,1));
                            set(this_FIP_signal.handles.data{i},...
                               'Visible','off',...
                               'XData',m_data{i}(:,2),...
                               'YData',m_data{i}(:,1));

                           % Ylabel
                           if use_current || same_y_label.raw
                               this_FIP_signal.handles.ylabel.String=...
                                   this_FIP_signal.settings.signal_units{this_FIP_signal.c_signal,1};
                           else
                               this_FIP_signal.handles.ylabel.String='mixed signals'; 
                           end
                        else
                            set(this_FIP_signal.handles.sig_CD{i},...
                                'Visible','off',...
                                'XData',m_sig_CD{i}(:,2),...
                                'YData',m_sig_CD{i}(:,1));
                            set(this_FIP_signal.handles.sig_405{i},...
                                'Visible','off',...
                                'XData',m_sig_405{i}(:,2),...
                                'YData',m_sig_405{i}(:,1));
                            set(this_FIP_signal.handles.data{i},...
                               'Visible',visibility,...
                               'XData',m_data{i}(:,2),...
                               'YData',m_data{i}(:,1));
                           % Ylabel
                           if use_current || same_y_label.data
                               this_FIP_signal.handles.ylabel.String=...
                                   this_FIP_signal.settings.signal_units{this_FIP_signal.c_signal,2};
                           else
                               this_FIP_signal.handles.ylabel.String='mixed signals'; 
                           end     
                        end      
                    end

                    % LogAI plot
                    % The reason for this awkward system rather than just
                    % making unnecesary plots invisible is that we don't even
                    % want the LogAI getter function to load al available raw
                    % data, because there could be huge amounth of data points
                    % there. This way LogAI only copies selected channels from
                    % Raw data. (See the getter function for clarification).
                    % However, if a figure is presented before logAI channels
                    % are selected all channels are plotted and we are not
                    % going to move those around either.
                    dims=size(logAI_plot);
                    j=1; % This is to skip deleted logAI plots
                    for i=2:dims(2) % first column is the timeline
                        j=j+1;
                        while ~this_FIP_signal.log_AI_used(j)
                            j=j+1;
                        end
                        this_FIP_signal.handles.logAI_plots{j}.XData=logAI_plot(:,1);
                        this_FIP_signal.handles.logAI_plots{j}.YData=logAI_plot(:,i);
                        % Only display them if not timestamps display
                        if this_FIP_signal.AI_plots
                            this_FIP_signal.handles.logAI_plots{j}.Visible='on';
                        else
                            this_FIP_signal.handles.logAI_plots{j}.Visible='off';
                        end
                    end
                    
                    % Update time stamps plots
                    if isfield(this_FIP_signal.handles,'time_stamp_plots') %there are timestamps
                        for i=1:length(this_FIP_signal.handles.time_stamp_plots)
                            this_FIP_signal.handles.time_stamp_plots{i}.XData=this_FIP_signal.timestamps{i};
                            if this_FIP_signal.AI_plots
                                this_FIP_signal.handles.time_stamp_plots{i}.Visible='off';
                            else
                                this_FIP_signal.handles.time_stamp_plots{i}.Visible='on';
                            end
                        end
                    end
                    
                    % Work on the x_axis, use window and c_time
                    if strcmp(this_FIP_signal.window,'all time')
                        x_limits=[0 max(this_FIP_signal.info.max_time)];
                    else
                        m_c_time=this_FIP_signal.c_time;
                        m_window=this_FIP_signal.window;
                        x_limits=[m_c_time-0.5*m_window m_c_time+0.5*m_window];
                        this_FIP_signal.handles.time_bar.Value=m_c_time;
                    end
                    this_FIP_signal.handles.FIP_plot.XLim=x_limits;
                    this_FIP_signal.handles.logAI_plot.XLim=x_limits;

                    % Indicator text
                    this_FIP_signal.handles.text_1.String=this_FIP_signal.info.names{this_FIP_signal.c_signal};
                    
                    % Button strings
                    % Raw or norm. data button
                    if this_FIP_signal.raw_signal
                        this_FIP_signal.handles.norm_data_button.String='norm. data';
                    else
                        this_FIP_signal.handles.norm_data_button.String='raw signal';
                    end
                    
                    % Time stamps or log AI button
                    if this_FIP_signal.AI_plots
                        this_FIP_signal.handles.time_stamp_button.Value=0;
                        this_FIP_signal.handles.time_stamp_button.String='Time stamps';
                    else
                        this_FIP_signal.handles.time_stamp_button.Value=1;
                        this_FIP_signal.handles.time_stamp_button.String='LogAI data';
                    end
                    

                case 'new_time' % it's all the way up there
                    % Work on the x_axis
                        if strcmp(this_FIP_signal.window,'all time')
                            x_limits=[0 max(this_FIP_signal.info.max_time)];
                        else
                            m_c_time=this_FIP_signal.c_time;
                            m_window=this_FIP_signal.window;
                            x_limits=[m_c_time-0.5*m_window m_c_time+0.5*m_window];
                            this_FIP_signal.handles.time_bar.Value=m_c_time;
                        end
                        this_FIP_signal.handles.FIP_plot.XLim=x_limits;
                        this_FIP_signal.handles.logAI_plot.XLim=x_limits;
                otherwise
                    warning('Not clear about plot update type')
            end
            
            % Either way, the c_time indicator line is always updated
            temp=this_FIP_signal.handles.FIP_plot.YLim;
            m_c_time=this_FIP_signal.c_time;
            set(this_FIP_signal.handles.c_time_line,...
                'XData',[m_c_time m_c_time],...
                'YData',[temp(1) temp(2)]);
        end
        
        function close_req(this_FIP_signal, scr, ev)
            % deals with closing of the object
            this_FIP_signal.figure_open=false;
            delete(this_FIP_signal.handles.main_figure)
        end
        
        function key_press(this_FIP_signal, scr, ev)
            % Deals with keypresses in all windows
            
            % General info about key pressed (used by other functions)
            this_FIP_signal.key_down=true;
            this_FIP_signal.key_pressed=ev.Key;
            
            % Deal with keypress actions
            switch scr.Tag
                case 'main_figure'
                    switch ev.Key
                        case 'rightarrow'
                            if this_FIP_signal.c_signal<this_FIP_signal.np_signals
                                this_FIP_signal.c_signal=this_FIP_signal.c_signal+1;
                                this_FIP_signal.update_plots;
                            end
                        case 'leftarrow'
                            if this_FIP_signal.c_signal>1
                                this_FIP_signal.c_signal=this_FIP_signal.c_signal-1;
                                this_FIP_signal.update_plots;
                            end
                        case 'uparrow'
                            this_FIP_signal.handles.norm_data_button.Value=...
                                  ~this_FIP_signal.handles.norm_data_button.Value;
                            this_FIP_signal.raw_signal=~this_FIP_signal.raw_signal;
                            this_FIP_signal.update_plots;
                        case 'downarrow'
                            this_FIP_signal.handles.norm_data_button.Value=...
                                  ~this_FIP_signal.handles.norm_data_button.Value;
                            this_FIP_signal.raw_signal=~this_FIP_signal.raw_signal;
                            this_FIP_signal.update_plots;
                        otherwise
                            % Unrecognized key, but that is fine. Other
                            % functions might still use it.
                    end
                otherwise
                    warning('unknown objects produces keypresses')
            end
        end
        
        function key_release(this_FIP_signal, scr, ev)
            % deals with key release
            this_FIP_signal.key_down=false;
            this_FIP_signal.key_pressed='empty';
  
        end
        
        function make_context(this_FIP_signal)
            % Function is responsible for making all the context menus
            
            % FIP_plot context menu
            this_FIP_signal.handles.drop_down.FIP_plot=uicontextmenu;
            uimenu(this_FIP_signal.handles.drop_down.FIP_plot,...
                'Tag','FIP_plot',...
                'Label','Reset zoom',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu(this_FIP_signal.handles.drop_down.FIP_plot,...
                'Tag','FIP_plot',...
                'Label','Restore raw data',...
                'Callback',@this_FIP_signal.context_menu)
            
            % Data context menu
            this_FIP_signal.handles.drop_down.data=uicontextmenu;
            temp_menu=uimenu(this_FIP_signal.handles.drop_down.data,...
                'Label','smooth_Trace');
            uimenu('Parent',temp_menu,'Tag','data',...
                'Label','smooth 1 sec',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu('Parent',temp_menu,'Tag','data',...
                'Label','custom',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu('Parent',temp_menu,'Tag','data',...
                'Label','no smooth',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu(this_FIP_signal.handles.drop_down.data,...
                'Label','export data',...
                'tag','data',...
                'Callback',@this_FIP_signal.context_menu);
            temp_menu=uimenu(this_FIP_signal.handles.drop_down.data,...
                'Label','data type');
            uimenu('Parent',temp_menu,'Tag','data',...
                'Label','deltaF/F',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu('Parent',temp_menu,'Tag','data',...
                'Label','F',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu('Parent',temp_menu,'Tag','data',...
                'Label','Z-score',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu(this_FIP_signal.handles.drop_down.data,...
                'Label','peri-event plot',...
                'Tag','data',...
                'Callback',@this_FIP_signal.context_menu);
            
            % 405nm context menu
            this_FIP_signal.handles.drop_down.sig_405=uicontextmenu;
            temp_menu=uimenu(this_FIP_signal.handles.drop_down.sig_405,...
                'Label','smooth_trace');
            uimenu('Parent',temp_menu,'Tag','sig_405',...
                'Label','smooth 1 sec',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu('Parent',temp_menu,'Tag','sig_405',...
                'Label','custom',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu('Parent',temp_menu,'Tag','sig_405',...
                'Label','no smooth',...
                'Callback',@this_FIP_signal.context_menu)
            temp_menu=uimenu(this_FIP_signal.handles.drop_down.sig_405,...
                'Label','fit method');
            uimenu('Parent',temp_menu,'Tag','sig_405',...
                'Label','polyfit_2 (standard)',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu('Parent',temp_menu,'Tag','sig_405',...
                'Label','polyfit_1',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu('Parent',temp_menu,'Tag','sig_405',...
                'Label','fit means',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu('Parent',temp_menu,'Tag','sig_405',...
                'Label','no fit',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu('Parent',temp_menu,'Tag','sig_405',...
                'Label','sliding window',...
                'Callback',@this_FIP_signal.context_menu)
            
            % Drop-down menus for logAI plot
            this_FIP_signal.handles.drop_down.logAI=uicontextmenu;
            temp_menu=uimenu(this_FIP_signal.handles.drop_down.logAI,...
                'Label','use for time allignment');
            uimenu('Parent',temp_menu,'Label','first peak is 5sec',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu('Parent',temp_menu,'Label','first peak is 15sec',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu('Parent',temp_menu,'Label','first peak is Xsec',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu(this_FIP_signal.handles.drop_down.logAI,...
                'Label','derive timestamps',...
                'Tag','log_AI',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu(this_FIP_signal.handles.drop_down.logAI,...
                'Label','delete',...
                'Tag','log_AI',...
                'Callback',@this_FIP_signal.context_menu)
            
            % Drop-down menu for timestamps
            this_FIP_signal.handles.drop_down.time_stamps=uicontextmenu;
            uimenu(this_FIP_signal.handles.drop_down.time_stamps,...
                'Label','derive timestamps',...
                'Tag','times_stamps',...
                'Callback',@this_FIP_signal.context_menu)
            uimenu(this_FIP_signal.handles.drop_down.time_stamps,...
                'Label','peri-event plot',...
                'Tag','times_stamps',...
                'Callback',@this_FIP_signal.context_menu) 
        end
        
        function context_menu(this_FIP_signal, scr, ev)
            % Deals with all context menus. I realize it's quite a long
            % list, but I've found it usefull to put all context menu
            % options in one big function. They are quite organised and
            % it's possible to use cmd-F (ctr-F) to search any context menu
            % label.
            
            % Note: when adding new context menu options. Add the option
            % itself to the correct menu in the function make_context.
            % Reffer to this function in the callback property. It migh be
            % usefull to know that the property this_FIP_signal.r_mouse_scr
            % contains the source of the last right mouse button click.

            % If the context menu option is simple and not used for
            % multiple options I try to put all the code in here rather
            % then making a seperate function.
            
            % Scroll to possible actions
            switch scr.Label
                case 'smooth 1 sec'
                    scr=this_FIP_signal.r_mouse_scr;
                    switch scr.Tag(1:5)
                        case 'sig_4'
                            this_FIP_signal.settings.smooth_405=this_FIP_signal.info.framerate;
                        case 'sig_C'
                            this_FIP_signal.settings.smooth_CD=this_FIP_signal.info.framerate;
                        case 'data.'
                            this_FIP_signal.settings.smooth_data=this_FIP_signal.info.framerate;
                        otherwise
                            warning('Unclear which signal should be smoothed')
                    end
                case 'custom'
                    scr=this_FIP_signal.r_mouse_scr;
                    smooth_value=inputdlg('# smooth datapoints');
                     switch scr.Tag(1:5)
                        case 'sig_4'
                            this_FIP_signal.settings.smooth_405=str2num(smooth_value{1});
                        case 'sig_C'
                            this_FIP_signal.settings.smooth_CD=str2num(smooth_value{1});
                        case 'data.'
                             this_FIP_signal.settings.smooth_data=str2num(smooth_value{1});
                        otherwise
                            warning('Unclear which signal should be smoothed')
                     end
                case 'no smooth'
                    smooth_value=1;
                    scr=this_FIP_signal.r_mouse_scr;
                     switch scr.Tag(1:5)
                        case 'sig_4'
                            this_FIP_signal.settings.smooth_405=smooth_value;
                        case 'sig_C'
                            this_FIP_signal.settings.smooth_CD=smooth_value;
                        case 'data.'
                             this_FIP_signal.settings.smooth_data=smooth_value;
                        otherwise
                            warning('Unclear which signal should be un-smoothed')
                     end
                case 'polyfit_2 (standard)'
                    this_FIP_signal.settings.fit_405='polyfit_2';
                case 'polyfit_1'
                    this_FIP_signal.settings.fit_405='polyfit_1';
                case 'fit means'
                    this_FIP_signal.settings.fit_405='fit_mean';
                case 'no fit'
                    this_FIP_signal.settings.fit_405='unfit';
                case 'sliding window'
                    this_FIP_signal.settings.fit_405='sliding_window';
                case 'export data'
                    variable_name=inputdlg('Data name:');
                    assignin('base',variable_name{1},...
                        this_FIP_signal.data{this_FIP_signal.c_signal});
                case 'Reset zoom'
                    this_FIP_signal.window='all time';
                    this_FIP_signal.update_plots;
                case 'Restore raw data'
                    % Restoring raw data (mostly settings)
                    m_settings=this_FIP_signal.settings;
                    m_settings.time_offset=0;
                    m_settings.smooth_405=1;
                    m_settings.smooth_CD=1;
                    m_settings.smooth_data=1;
                    m_settings.fit_405='unfit';
                    for i=1:this_FIP_signal.np_signals
                        m_settings.signal_units{i,1}='F';
                        m_settings.signal_units{i,2}='deltaF/F';
                    end
                    this_FIP_signal.settings=m_settings;
                    % Other:
                    this_FIP_signal.log_AI_used=...
                        true(length( this_FIP_signal.log_AI_used),1);
                    % Update
                    this_FIP_signal.update_plots;
                case 'Z-score'
                    scr=this_FIP_signal.r_mouse_scr;
                    if strcmp(scr.Tag(1:4),'data')
                        s_type=2; % Data
                    else % This needs work
                        s_type=1; % CD or 405nm signal
                    end
                    index=find('.'==scr.Tag,1)+1; %after the point is the signal #
                    s_np=str2num(scr.Tag(index:end));
                    this_FIP_signal.settings.signal_units{s_np,s_type}='Z-score';
                case 'deltaF/F'
                    scr=this_FIP_signal.r_mouse_scr;
                    if strcmp(scr.Tag(1:4),'data')
                        s_type=2; % Data
                    else % This needs work
                        s_type=1; % CD or 405nm signal
                    end
                    index=find('.'==scr.Tag,1)+1; %after the point is the signal #
                    s_np=str2num(scr.Tag(index:end));
                    this_FIP_signal.settings.signal_units{s_np,s_type}='deltaF/F';
                case 'F'
                    scr=this_FIP_signal.r_mouse_scr;
                    if strcmp(scr.Tag(1:4),'data')
                        s_type=2; % Data
                    else % This needs work
                        s_type=1; % CD or 405nm signal
                    end
                    index=find('.'==scr.Tag,1)+1; %after the point is the signal #
                    s_np=str2num(scr.Tag(index:end));
                    this_FIP_signal.settings.signal_units{s_np,s_type}='F';
                case 'delete'
                    % delete src of R mouse click
                    scr=this_FIP_signal.r_mouse_scr;
                    if strcmp(scr.Tag(1:3),'AI ')
                        % It's a logAI trace
                        AI_np=str2num(scr.Tag(4:end));
                        this_FIP_signal.log_AI_used(AI_np)=false;
                    end 
                    delete(scr) %deletes the actual plot, not super ellegant
                case 'derive timestamps'
                    % Derive timestamps from input
                    tag=this_FIP_signal.r_mouse_scr.Tag;
                    if strcmp(tag(1:4),'stam')
                        % Deriving stamps from stamps
                        nr=str2num(tag(8:end));
                        input=inputdlg('minimum stamp interval: ');
                        input=str2num(input{1});
                        stamps_old=this_FIP_signal.timestamps{nr};
                        j=0;
                        j=1;
                        stamps(j)=stamps_old(1);
                        for i=2:length(stamps_old)
                            if stamps_old(i)-stamps_old(i-1)>input
                                j=j+1;
                                stamps(j)=stamps_old(i);
                            end
                        end
                        name=[this_FIP_signal.timestamps_names{nr} '_meta'];
                    else
                        stamps=this_FIP_signal.derive_stamps(this_FIP_signal.r_mouse_scr);
                        name=[tag '_derived'];
                    end
                    % Store them in the FIP_signal object and name them
                    this_FIP_signal.import_timestamps(stamps,name);
                case 'peri-event plot'
                    % Make peri-event data based on scr
                    tag=this_FIP_signal.r_mouse_scr.Tag;
                    if strcmp(tag(1:4),'stam')
                        nr=str2num(tag(8:end));
                        stamps=this_FIP_signal.timestamps{nr};
                        qstring=['Usign the signal from '...
                        this_FIP_signal.info.names{this_FIP_signal.c_signal}];
                        
                        input=questdlg(qstring,'Select data',...
                            'Calcium Dependend',...
                            '405nm signal',...
                            'normalized signal',...
                            'normalized signal');

                        switch input
                            case 'Calcium Dependend'
                                m_data=this_FIP_signal.sig_CD{this_FIP_signal.c_signal};
                            case '405nm signal'
                                m_data=this_FIP_signal.sig_405{this_FIP_signal.c_signal};
                            case 'normalized signal'
                                m_data=this_FIP_signal.data{this_FIP_signal.c_signal};
                        end
                    elseif strcmp(tag(1:4),'data') %clicked on data
                        input=listdlg('PromptString','Select time stamps:',...
                            'SelectionMode','single',...
                            'ListString',this_FIP_signal.timestamps_names);
                        stamps=this_FIP_signal.timestamps{input};
                        nr=str2num(tag(6:end));
                        m_data=this_FIP_signal.data{nr};
                    end
                    
                    
                   
                    window=15; %Will work on changing this later
                    this_FIP_signal.peri_event_plot(m_data,stamps,window);
                case 'first peak is 5sec'
                    % Find time of first peak
                    stamps=this_FIP_signal.derive_stamps(this_FIP_signal.r_mouse_scr);
                    time_offset=stamps(1)-5;
                    disp(['Time offset is ' num2str(time_offset) ' sec.'])
                    this_FIP_signal.settings.time_offset=time_offset;
                    this_FIP_signal.update_plots;
                    % Check if there are allready timestamps stored
                    if isfield(this_FIP_signal.handles,'time_stamp_plots')
                        warning('There all allready timestamps, those will NOT be offset.')
                    end
                case 'first peak is 15sec'
                    % Find time of first peak
                    stamps=this_FIP_signal.derive_stamps(this_FIP_signal.r_mouse_scr);
                    time_offset=stamps(1)-15;
                    disp(['Time offset is ' num2str(time_offset) ' sec.'])
                    this_FIP_signal.settings.time_offset=time_offset;
                    this_FIP_signal.update_plots;
                    % Check if there are allready timestamps stored
                    if isfield(this_FIP_signal.handles,'time_stamp_plots')
                        warning('There all allready timestamps, those will NOT be offset.')
                    end
                otherwise
                    warning('Menu option not yet available.')
            end
        end
        
        function mouse_down(this_FIP_signal, scr, ev)
            % Deals with mouse clicks

            % Currently only source is FIP plot, so not dealing with scr
            if ev.Button==1
                % clicked on line, not on plot -> set Parents as source
                if strcmp(scr.Type,'line')
                    scr=scr.Parent;
                end
                
                m_mouse_action='none'; %please leave it a 4-char.
                if this_FIP_signal.key_down
                    if strcmp(this_FIP_signal.key_pressed,'shift')
                        m_mouse_action='drag';
                    end
                end 
                
                % store initial click location
                punter=get(scr,'CurrentPoint');
                this_FIP_signal.mouse_click=true;
                this_FIP_signal.mouse_start=punter(1,1:2);
                
                if m_mouse_action=='drag'
                    % Mouse button drag/scroll
                    punter=get(scr,'CurrentPoint');
                    this_FIP_signal.mouse_click=true;
                    this_FIP_signal.mouse_start=punter(1,1:2);
                    this_FIP_signal.mouse_action.name='time_scroll';
                    this_FIP_signal.mouse_action.start_time=...
                        this_FIP_signal.c_time;
                    this_FIP_signal.mouse_action.total_time_dif=0;
                    % We don't want to evaluate every little mouse
                    % movement, so we implement a counter.
                    this_FIP_signal.mouse_action.counter=0;
                    this_FIP_signal.mouse_action.scr=scr; 
                else
                    % Horizontal zoom is the only other mouse action rn
                    this_FIP_signal.mouse_action.name='hor_zoom';
                    this_FIP_signal.mouse_action.drawing=...
                        line([punter(1,1) punter(1,1)],[punter(1,2) punter(1,2)],...
                        'Color',[0 0 0]);
                    this_FIP_signal.mouse_action.scr=scr;
                end
                
                % Set the mouse motion Fcn for scr parent (axes don't
                % have motionFcn or button up Fcn, but figures do.
                scr.Parent.WindowButtonMotionFcn=@this_FIP_signal.mouse_motion;
                scr.Parent.WindowButtonUpFcn=@this_FIP_signal.mouse_release;
                
            elseif ev.Button==3
                % Store the source of the R click for the context menu
                this_FIP_signal.r_mouse_scr=scr;
            end
        end
        
        function mouse_motion(this_FIP_signal, scr, ev)
            % Deals with mouse motion
            
            % Find out the type of action
            switch this_FIP_signal.mouse_action.name
                case 'hor_zoom'
                    % Currently the only action
                    scr=this_FIP_signal.mouse_action.scr;
                    punter=scr.CurrentPoint;
                    this_FIP_signal.mouse_action.drawing.XData(2)=punter(1,1);
                case 'time_scroll'
                    this_FIP_signal.mouse_action.counter=...
                        this_FIP_signal.mouse_action.counter+1;
                    if this_FIP_signal.mouse_action.counter==1 % Will update every n mouse movements
                        this_FIP_signal.mouse_action.counter=0;
                        scr=this_FIP_signal.mouse_action.scr;
                        punter=scr.CurrentPoint;
                        time_diff=punter(1,1)-this_FIP_signal.mouse_start(1,1);
                        this_FIP_signal.c_time=this_FIP_signal.c_time-time_diff;
                        this_FIP_signal.mouse_action.total_time_dif=...
                            this_FIP_signal.mouse_action.total_time_dif+time_diff;
                        this_FIP_signal.update_plots('type','new_time');
                    end
                otherwise
                    warning('mouse error han_114')
            end
        end
        
        function mouse_release(this_FIP_signal, scr, ev)
            % Deals with mouse release
            
            % Get the mouse click source
            scr=this_FIP_signal.mouse_action.scr;
            
            % Find out what action was performed
            switch this_FIP_signal.mouse_action.name
                case 'hor_zoom'
                    % Set xaxis to new limits
                    punter=scr.CurrentPoint;
                    if abs(punter(1,1)-this_FIP_signal.mouse_start(1))>1
                        if punter(1,1)<this_FIP_signal.mouse_start(1)
                            new_XLim=[punter(1,1),...
                                this_FIP_signal.mouse_start(1)];
                        else
                            new_XLim=[this_FIP_signal.mouse_start(1),...
                                punter(1,1)]; % Have to deal with negative lines
                        end
                        this_FIP_signal.c_time=mean(new_XLim);
                        this_FIP_signal.window=new_XLim(2)-new_XLim(1);
                        this_FIP_signal.update_plots;
                    end
                    delete(this_FIP_signal.mouse_action.drawing);
                case 'time_scroll'
                    % well stop scrolling 
                    all_good=true;
                otherwise
                    warning('mouse error han_115')
            end
            
            % Delete mouse action and Fcn in figure
            scr.Parent.WindowButtonMotionFcn='';
            scr.Parent.WindowButtonUpFcn='';
            this_FIP_signal.mouse_click=false;
            this_FIP_signal.mouse_action.name='none';
            % Leaving the original click in mouse_start
        end
        
        function UI_used(this_FIP_signal, scr, ev)
            % This function deals with all UI ellement interaction
            
            %assignin('base','testera',scr)
            button=scr.Tag(1:6); %because put random stuff behind
            
            switch button
                case 'signal'
                    index=find('_'==scr.Tag,1)+1;
                    i=str2num(scr.Tag(index:end));
                    this_FIP_signal.fibers_shown(i)=logical(scr.Value);
                    this_FIP_signal.update_plots;
                case 'data_s'
                    if this_FIP_signal.raw_signal
                        scr.String='raw signal';
                    else
                        scr.String='norm. data';
                    end
                        
                    this_FIP_signal.raw_signal=~this_FIP_signal.raw_signal;
                    this_FIP_signal.update_plots;
                case 'time_s'
                    m_c_time=scr.Value;
                    this_FIP_signal.c_time=m_c_time;
                    this_FIP_signal.update_plots('type','new_time');
                case 'stamps'
                    this_FIP_signal.AI_plots=~this_FIP_signal.AI_plots;
                    this_FIP_signal.update_plots;
                otherwise
                    disp(button)
                    disp('Function currently not suported')
            end  
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    events
       state_change         % fires at every change (triggerd by update_plots)
    end
    
end

