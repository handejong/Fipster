classdef FIP_acquisition < handle
    %FIP_ACQUISITION Aquires multi-channel fiber photometry data
    %   FIP_acquisition objects are used to acquire multi-channel FIP data
    %   using a camera. The camera takes snapshots of a fiber bundle at a
    %   certain frequency. Individual fibers are then segmented online.
    %
    %   The aquisition method is based on FIPGUI by the Deisseroth lab. See
    %   Kim et al., Nat. Methods 2016 for rationale and Github link.
    %
    %   Currently supported DAQ:
    %       - National Instruments
    %
    %   Currently supported cameras:
    %       - Photometrics Prime (note: Matlab adapter in B�ta)
    %       - Macvideo facetime camera
    %
    %   FIP_acquisition is part of FIPster. FIPster is made by Johannes de
    %   Jong, j.w.dejong@berkeley.edu
    
    
    properties
        filename        % How the saved datafiles wil be named
        cam_settings    % Struct with camera settings and video source object
        daq_settings    % Struct with DAQ info and session
        aq_settings     % Struct with acquisition settings
        plt_settings    % Struct with plot settings
        handles         % Mostly plot handles
        notes           % Cell struct with notes taken during recordings
    end
    
    methods
        function obj = FIP_acquisition(varargin)
            % FIP_ACQUISITION Construct an instance of this class
            %   Will connect to harware and populate settings
            
            % Input argument handeling
                % TODO: including GUI
                
            % Default variables
            rate = 20; % Hz (in total!)
            sample_rate_factor = 10; % how much faster DAQ samples than camera
            exposure_gap = 10; % ms (Used for reading and cleaning of sensor)
            lookback = 10; % s plotting signal history
            backup_cycle = 0; % Block of datapoints writen to csv during acquisition (0 = no backup)
              
            % Find all camera options
            obj.getCameraHardware();
             
            % User will be able to set camera, but pick one for now:
            obj.cam_settings.in_use=1;
            
            % Setup DAQ
            try
                obj.daq_settings.devices = daq.getDevices();
                obj.daq_settings.in_use = 1;
            catch
                disp('No data aquisition toolbox installed, or unable to connect to DAQ.')
                obj.daq_settings.devices = 'none';
                obj.daq_settings.in_use = 'none';
            end
            
            % Populate DAQ settings
            obj.daq_settings.fs = rate * sample_rate_factor;
            obj.daq_settings.SRF = sample_rate_factor;
            obj.daq_settings.ports.camera = 0;
            obj.daq_settings.ports.sig = 2;
            obj.daq_settings.ports.ref = 1;
            
            % Populate acquisition settings
            obj.aq_settings.calibrated = false;
            obj.aq_settings.rate = rate;
            obj.aq_settings.exposure_gap = exposure_gap;
            obj.aq_settings.channels = 2;
            obj.aq_settings.fibers = 1;
            obj.aq_settings.acquire_AI = true; % Acquire AI signals as well
            obj.aq_settings.backup_cycle = backup_cycle;
            
            % Populate plot settings
            obj.plt_settings.lookback=  lookback;
            obj.plt_settings.type = '.';
            obj.plt_settings.smooth = 1;
            obj.plt_settings.live_feed = false; 
            
            % Store empty notes
            obj.notes={'These notes are saved after recording.'};
        end
        
        function calibrate(obj)
            %calibraTE will calibrate hardware
            %   Presents snapshot that should be used to select fibers,
            %   then runs program to set LED's to same signal strength.
            
            % Default variables
            n_frames = 4; %number for frames used for calibration
            
            % Setup the selected camera
            obj.camera_setup;
            
            % Set the ROI to the whole screen
            res=obj.cam_settings.vid.VideoResolution;
            obj.cam_settings.vid.ROIPosition=[0 0 res];
            
            
            % If neccesary setup daq
            if strcmp(obj.cam_settings.trigger_type, 'daq')
                obj.daq_setup;
            end
            
            % Where we will store calibration frames
            frames = zeros(res(2), res(1), n_frames);
            
            if ~strcmp(obj.cam_settings.trigger_type,'daq')
                start(obj.cam_settings.vid);
                pause(1)
                i=0;
                while i < n_frames
                    i = i + 1;
                    trigger(obj.cam_settings.vid)
                    while(islogging(obj.cam_settings.vid)); end
                    frames(:,:,i) = getdata(obj.cam_settings.vid, 1, 'uint16');
                end
                stop(obj.cam_settings.vid);
            else
                start(obj.cam_settings.vid);
                startBackground(obj.daq_settings.session);
                i=0;
                while i < n_frames
                    i = i + 1;
                    frames(:,:,i) = getdata(obj.cam_settings.vid, 1, 'uint16');
                end
                stop(obj.cam_settings.vid);
                stop(obj.daq_settings.session);   
            end
            
            % Make the four frames into one calibframe
            calibframe = max(frames, [], 3);
            
            % Fiber ROI GUI, update acquisition settings.
            calibOut = calibrationgui(calibframe);
            obj.aq_settings.masks = calibOut.masks;
            obj.aq_settings.calibImg = calibOut.figImg;
            obj.aq_settings.labels = calibOut.labels;
            obj.aq_settings.fibers=length(calibOut.labels);
            
            % Set camera ROI based on where the fibers are
            ROI=[0 0 0 0];
            super_mask=sum(obj.aq_settings.masks,3);
            super_mask(super_mask~=0)=1;
            X=sum(super_mask,1);
            Y=sum(super_mask,2);
            i=0; temp=0;
            while temp==0; i=i+1; temp=X(i); end; ROI(1)=i-1;
            i=length(X)+1; temp=0;
            while temp==0; i=i-1; temp=X(i); end; ROI(3)=(i+1)-ROI(1);
            i=0; temp=0;
            while temp==0; i=i+1; temp=Y(i); end; ROI(2)=i-1;
            i=length(Y)+1; temp=0;
            while temp==0; i=i-1; temp=Y(i); end; ROI(4)=(i+1)-ROI(2);
            obj.cam_settings.vid.ROIPosition=ROI;
            
            % Update the masks themselves
            obj.aq_settings.masks=logical(obj.aq_settings.masks(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1,:));
            
            % Update the calibrated for setting
            obj.cam_settings.calibrated_for=obj.cam_settings.in_use;
            
            % Set aqcuisition boolean to 'is-calibrated'
            obj.aq_settings.calibrated = true;
        end
        
        function start_acquisition(obj)
            %START_ACQUISITION does the actual data acquisition
            
            % First check if properly calibrated
            
            % %%% TO DO! %%%
            
            % Check the filename, and increment if neccesary
            % This could be a separate function
            if isempty(obj.filename)
                temp=inputdlg('Please give a valid filename: ');
                obj.filename=temp{1};
            end
            while exist([obj.filename '.mat'],'file')==2
                obj.filename=[obj.filename '_copy'];
            end
            
            
            % Acquisition variables
            n_frame = 0;
            data = zeros(obj.aq_settings.fibers, 1);
            signal = zeros(100, 2, obj.aq_settings.fibers);
            ref = zeros(size(signal));
            lookback = (obj.plt_settings.lookback * obj.aq_settings.rate) / obj.aq_settings.channels;
            
           
            % Plot
            main_figure = figure;
            for i = 1:obj.aq_settings.fibers
                subplot(obj.aq_settings.fibers,1,i)
                yyaxis left
                signal_plots{i} = plot(signal(:,2,i), signal(:,1,i), obj.plt_settings.type);
                hold on
                axis tight
                ylabel('Signal (F)')
                if obj.aq_settings.channels == 2
                    yyaxis right
                    ref_plots{i} = plot(ref(:,2,i), ref(:,2,i), obj.plt_settings.type);
                    ylabel('Refference (F)')
                end
                xlabel('time (s)')
            end
            
            
            % Make a figure that shows a live feed of the fibers
            if obj.plt_settings.live_feed
                figure
                live_feed = imagesc(max(double(obj.aq_settings.masks),[],3));
                caxis([200 250])
            end
            
            
            % Make a figure for the stop button and notepad
            stop_figure = figure('Menubar','none',...
                'CloseRequestFcn', @obj.uncloseable);
            stop_figure.Position(3:4) = [300, 400];
            
            
            % Stop button
            stop_button = uicontrol(stop_figure,...
                'Units','normalized',...
                'Position',[0.1 0.85 0.8 0.1],...
                'Style','togglebutton',...
                'String','Stop acquisition',...
                'HorizontalAlignment','center');
            
            % Notepad
            uicontrol(stop_figure,...
                'Style','text',...
                'String','Add notes below: ',...
                'Units','normalize',...
                'Position',[0.1 0.73 0.8 0.1],...
                'HorizontalAlignment','left');
            
            notes_input = uicontrol(stop_figure,...
                'Style','edit',...
                'Units','normalized',...
                'Position',[0.1 0.1 0.8 0.66],...
                'Max', 100,...
                'HorizontalAlignment','left',...
                'String',obj.notes);
            
            
            % Start the camera and daq depending on the camera type
            obj.camera_setup();
            start(obj.cam_settings.vid);
            if strcmp(obj.cam_settings.trigger_type, 'daq')
                obj.daq_setup();
                startBackground(obj.daq_settings.session)
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%% Acquisition loop %%%%%%%%%%%%%%%%%%%%%%%
            should_be_time = 0;
            start_time = tic;
            error_counter=0;
            analysis_time_timer = tic; % only the first time
            while(stop_button.Value == 0)
                
                
                % Frame number and expected time
                n_frame = n_frame+1;
                should_be_time = should_be_time + 1 / obj.aq_settings.rate;
                
                % Exponentially increase the size of the signal and ref
                % variables
                if n_frame>length(signal)
                    signal = [signal; zeros(size(signal))];
                    ref=[ref; zeros(size(ref))];
                end
                
                % Get analysis time (loop starte below acquisition)
                analysis_time = toc( analysis_time_timer );
                
                % Start acquisition timer
                acquisition_timer = tic;

                % Take a picture
                switch obj.cam_settings.trigger_type
                    case 'internal'
                        trigger(obj.cam_settings.vid)
                        while(islogging(obj.cam_settings.vid)); end
                        [img, time, metadata] = getdata(obj.cam_settings.vid,1);
                        
                    case 'daq'
                        try
                            [img, time, metadata] = getdata(obj.cam_settings.vid,1,'uint16');
                        catch
                            disp('No frames available in camera memory.')
                            stop_button.Value = 1;
                        end

                    otherwise
                        error('Incompatible trigger type.')
                end
                
                % Time acquisition
                acquisition_time = toc(acquisition_timer);
                
                % To calculate analysis time
                analysis_time_timer = tic;
                
                % Error handeling, time
                if toc(start_time) - should_be_time > 1
                    warning(['Acquisition is delayed by ' num2str(toc(start_time) - should_be_time) ' s'])
                    disp(['Acquisition time: ' num2str(acquisition_time) ' s'])
                    disp(['Analysis time: ' num2str(analysis_time) ' s'])
                end
                
                % Error handeling, frame number
                if metadata.FrameNumber~=n_frame
                    warning('Camera and acquisition seem to disagree on frame number')
                end
                
                % Segment the individual fibers out
                for j = 1:obj.aq_settings.fibers
                    data(j) = mean(img(obj.aq_settings.masks(:,:,j)));
                end
                
                % Store the data in the correct channel
                sr_pair = ceil(n_frame / obj.aq_settings.channels);
                if mod(n_frame,obj.aq_settings.channels) == 0
                    signal(sr_pair, 1, :) = data;
                    signal(sr_pair, 2, :) = toc(start_time);
                else
                    ref(sr_pair, 1, :) = data;
                    ref(sr_pair, 2, :) = toc(start_time);
                end
                
                % Update live feed (if requested)
                if obj.plt_settings.live_feed && mod(n_frame,obj.aq_settings.channels) == 0
                    live_feed.CData=img;
                end
                
                % Save data with raw timestamps in a .csv for backup
                % Note that the collums are time, all signals, time, all refs
                if mod(sr_pair, obj.aq_settings.backup_cycle) == 0 && mod(n_frame,obj.aq_settings.channels) == 0
                    index = obj.aq_settings.backup_cycle -1;
                    temp=[squeeze(signal(sr_pair-index:sr_pair, 2, 1)) squeeze(signal(sr_pair-index:sr_pair, 1,:)) squeeze(ref(sr_pair-index:sr_pair, 2, 1)) squeeze(ref(sr_pair-index:sr_pair, 1,:))];
                    dlmwrite([obj.filename '.csv'],...
                        temp,...
                        'delimiter',',','-append');
                end
                
                % Update plots
                if sr_pair > lookback
                    start_index = sr_pair-lookback;
                    for j = 1:obj.aq_settings.fibers
                        signal_plots{j}.XData = signal(start_index:sr_pair-1, 2, j);
                        signal_plots{j}.YData = smooth(signal(start_index:sr_pair-1, 1, j), obj.plt_settings.smooth);
                        if obj.aq_settings.channels == 2
                            ref_plots{j}.XData=ref(start_index:sr_pair-1, 2, j);
                            ref_plots{j}.YData=smooth(ref(start_index:sr_pair-1, 1, j), obj.plt_settings.smooth);
                        end
                    end
                else % plot signal from the beginning
                    for j = 1:obj.aq_settings.fibers
                        signal_plots{j}.XData=signal(1:sr_pair, 2, j);
                        signal_plots{j}.YData=signal(1:sr_pair, 1, j);
                        if obj.aq_settings.channels == 2
                            ref_plots{j}.XData=ref(1:sr_pair, 2, j);
                            ref_plots{j}.YData=ref(1:sr_pair, 1, j);
                        end
                    end
                end
                drawnow();
                
                % If not using a daq or pulse generator, pause to sync
                if strcmp(obj.cam_settings.trigger_type, 'internal')
                    if toc(start_time) < should_be_time
                        pause(should_be_time - toc(start_time));
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%% End of acquisition loop %%%%%%%%%%%%%%%%%%%
            
            % Check if we got all frames
            if obj.cam_settings.vid.FramesAvailable~=0
                warning(['Camera has ' num2str(obj.cam_settings.vid.FramesAvailable) ' frames available, this could be a problem.'])
            end    
            
            % Collect all the notes
            obj.notes = notes_input.String;
            
            % Delete  stop button
            close(stop_figure, 'force')
            
            % Stop the camera
            stop(obj.cam_settings.vid);
            
            % If neccesary, stop the DAQ
            if strcmp(obj.cam_settings.trigger_type,'daq')
                stop(obj.daq_settings.session)
            end
            
            % Delete unused signal
            sr_pair = floor(n_frame / obj.aq_settings.channels);
            signal = signal(1:sr_pair, :, :);
            ref = ref(1:sr_pair, :, :);
            
            % Update the timeline.
            %   Info: even altough we've stored timestamps in the signal
            %   and/or reference arrays, if we used a daq to trigger the
            %   camera and the LED's, these are probably more accurate and
            %   we should just trust those. It might be that the aquisition
            %   loop was lagging behind (frames were temporarily stored in
            %   camera memory). It is not possible that the aquisition loop
            %   went faster than the camera, because no frames would be
            %   available and the getdata function would fail.
            
            % The start times are 0.5*interval and 1.5 * interval, because the
            % first trigger is at 0, but the first readout is at
            % T = exposure.
            if strcmp(obj.cam_settings.trigger_type,'daq')
                interval=(1/obj.aq_settings.rate)*obj.aq_settings.channels;
                if obj.aq_settings.channels == 2; signal_start=0.75 * interval; end
                if obj.aq_settings.channels == 1; signal_start = 0.5 * interval; end
                ref_start=0.25 * interval;

                for i=1:obj.aq_settings.fibers
                    signal(:,2,i)=[signal_start:interval:length(signal)*interval];
                    ref(:,2,i)=[ref_start:interval:length(ref)*interval];
                end
                disp('Exact timestamps using daq timer.')
            end
            
            % NOTE: the backup data (saved as .csv if requested) contains
            % the original timestamps at which datapoints were collected
            % by the getdata function. This is the approximate time at
            % read-out, which is an integral of the signal between trigger
            % and read-out. The timestamps should be displaced by 1/4th of
            % the sampling rate (2 channels) or 1/2 the sampling rate (1
            % channel).
            
            % Save the data
            framerate = obj.aq_settings.rate;
            labels = obj.aq_settings.labels;
            notes = obj.notes;
            if obj.aq_settings.channels==1
                save(obj.filename,'signal','framerate','labels','notes','-v7.3')
            else
                save(obj.filename,'signal','ref','framerate','labels','notes','-v7.3')
            end
            imwrite(obj.aq_settings.calibImg.cdata, [obj.filename '_calibration.jpg'], 'JPEG') % calibration image

            % Plot all the data
            figure(main_figure)
            for i=1:obj.aq_settings.fibers
                signal_plots{i}.XData = signal(1:sr_pair, 2, i);
                signal_plots{i}.YData = signal(1:sr_pair, 1, i);
                ref_plots{i}.XData = ref(1:sr_pair, 2, i);
                ref_plots{i}.YData = ref(1:sr_pair, 1, i);
                
                % Plot smoothed data on top
                subplot(obj.aq_settings.fibers, 1, i)
                yyaxis left
                plot(signal(:,2,i), smooth(signal(:, 1, i), 10));
                yyaxis right
                plot(ref(:,2,i), smooth(ref(:, 1, i), 10));
            end  
            
            % Store the signal in the workspace
            assignin('base', 'new_signal', signal)
            assignin('base', 'new_ref', ref)
            
        end
        
        function snap = snap_frame(obj)
            %SNAP_FRAME snaps one frame using the current settings
            
            % setup camera
            obj.camera_setup;
            
            % snap picture
            snap = getsnapshot(obj.cam_settings.vid);
            
            % Present figure
            figure
            imagesc(snap);
            colorbar();
        end
       
    end
    
    methods (Access =private)
        
        function getCameraHardware(obj)
            %getCameraHardware list information about all attached cameras
            %   Returns four lists, with each row refering to an available camera
            %   device. This function is adapted from FIPGUI (Deisseroth lab).
            
            % Reset camera adaptors
            imaqreset();
            
            % Cycle to all available adaptors and check suport
            adaptorsOut = []; devicesOut = {}; formatsOut = []; IDsOut = [];
            i = 0;
            adaptors = imaqhwinfo();
            for adaptor = adaptors.InstalledAdaptors
                if obj.check_cam_support(adaptor)
                    devices = imaqhwinfo(adaptor{:});
                    for device = devices.DeviceInfo
                        for format = device.SupportedFormats
                            i = i + 1;
                            options(i).adaptors = adaptor{:};
                            devNameParts = strsplit(device.DeviceName, ',');
                            options(i).devices = devNameParts{1};
                            options(i).formats = format{:};
                            options(i).IDs = device.DeviceID;
                        end
                    end
                end
            end
            
            % Error if no supported cameras found
            if i==0
                error('MATLAB IMAQ detected no available camera devices to connect to. Fix this and restart MATLAB.');
            end
            
            % Import all adaptors and cameras into object
            obj.cam_settings.options=options;
            
        end
        
        function supported = check_cam_support(~, adaptor)
            %CHECK_CAM_SUPPORT tests if the current adaptor is supported
            
            switch adaptor{1}
                case 'macvideo'
                    supported = true;
                case 'pmimaq'
                    supported = true;
                otherwise
                    supported = false;
            end
            
            if supported
                disp([adaptor{1} ' is supported'])
            else
                disp([adaptor{1} ' is NOT supported'])
            end 
        end
        
        function camera_setup(obj)
            %CAMERA_SETUP checks the selected camera and sets it up
            
            % Figure out which camera the user wants to use
            camDeviceN = obj.cam_settings.in_use;
            
            
            % First check if allready setup
            if isfield(obj.cam_settings,'setup_for') && obj.cam_settings.setup_for==camDeviceN
                % don't have to remake object, but should still run other
                % settings because ROI and exposure might have changed.
                option=obj.cam_settings.options(camDeviceN);
                disp('No new video object')
            else
                option=obj.cam_settings.options(camDeviceN);
                vid = videoinput(option.adaptors, option.IDs, option.formats);
                src = getselectedsource(vid);
                obj.cam_settings.vid = vid;
                obj.cam_settings.src = src;
                disp(['Setting up ' option.devices])
            end
            
            
            % Figure out the correct settings for supported cameras
            switch option.adaptors
                case 'macvideo'
                    % Setup for Apple Facetime cam
                    obj.cam_settings.trigger_type='internal';
                    triggerconfig(obj.cam_settings.vid, 'Manual')
                    obj.cam_settings.vid.TriggerRepeat = Inf;
                    obj.cam_settings.vid.FramesPerTrigger=1;
                    obj.cam_settings.vid.ReturnedColorspace = 'grayscale';
                    
                case 'pmimaq'
                    % Setup for Photometrics Prime;
                    obj.cam_settings.vid.FramesPerTrigger = 1; 
                    obj.cam_settings.vid.TriggerRepeat = Inf;
                    
                    % Camera settings
                    obj.cam_settings.src.TriggerMode='Edge Trigger';
                    obj.cam_settings.src.AutoContrast = 'OFF';
                    obj.cam_settings.src.PP0ENABLED = 'NO';
                    obj.cam_settings.src.PP1ENABLED = 'NO';
                    obj.cam_settings.src.PP2ENABLED = 'NO';
                    obj.cam_settings.src.PP3ENABLED = 'NO';
                    obj.cam_settings.src.PP4ENABLED = 'NO';
                    obj.cam_settings.src.ClearCycles = 1;
                    obj.cam_settings.src.ClearMode = 'Pre-Exposure';
                    
                    % Set the ROI
                    %obj.cam_settings.vid.ROIPosition = [0 0 obj.cam_settings.vid.VideoResolution];
                    
                    % Set the exposure
                    obj.cam_settings.src.Exposure = (1000/obj.aq_settings.rate) - obj.aq_settings.exposure_gap;
                    
                    % set trigger type
                    obj.cam_settings.trigger_type = 'daq'; % use daq to trigger
                    
                    disp('Setup for Prime.')
                    
                otherwise
                    error([option.adaptors ' currently not suported. See cam_setup method.'])
            end
            
            % Update what for which camera the object is currently set
            obj.cam_settings.setup_for=camDeviceN;
        end
        
        function daq_setup(obj)
            %DAQ_SETUP sets upt the daq according to the settings
            
            % Get device from object
            device = obj.daq_settings.devices(obj.daq_settings.in_use);
            
            % Info about the ports
            rate=obj.aq_settings.rate;
            cam_port=obj.daq_settings.ports.camera;
            sig_port=obj.daq_settings.ports.sig;
            ref_port=obj.daq_settings.ports.ref;
            
            % Make a session
            s = daq.createSession('ni');
            s.Rate = obj.daq_settings.fs;
            s.IsContinuous = true;
            
            % Setting up the camera triggers
            camCh = s.addCounterOutputChannel(device.ID, cam_port, 'PulseGeneration');
            camCh.Frequency = rate;
            camCh.InitialDelay = 0;
            camCh.DutyCycle = 0.1;
            disp(['Camera should be connected to ' camCh.Terminal]);
            
            % Setting up the LED triggers (1 or 2 LEDs)
            switch obj.aq_settings.channels
                case 1
                    sigCh = s.addCounterOutputChannel(device.ID, sig_port, 'PulseGeneration');
                    sigCh.Frequency = rate;
                    sigCh.InitialDelay = 0;
                    sigCh.DutyCycle = 0.95;
                    disp(['Signal LED should be connected to ' sigCh.Terminal]);
                    
                case 2
                    refCh = s.addCounterOutputChannel(device.ID, ref_port, 'PulseGeneration');
                    refCh.Frequency = rate / 2;
                    refCh.InitialDelay = 1 / rate * 0.05;
                    refCh.DutyCycle  = 0.45;
                    disp(['Reference LED should be connected to ' refCh.Terminal]);

                    sigCh = s.addCounterOutputChannel(device.ID, sig_port, 'PulseGeneration');
                    sigCh.Frequency = rate / 2;
                    sigCh.InitialDelay = 1 / rate * 1.05;
                    sigCh.DutyCycle = 0.45;
                    disp(['Signal LED should be connected to ' sigCh.Terminal]);
                
                otherwise
                    error('Currently only setup for 1 and 2 channel recordings.')
            end
            
            % Enabling analog input (AI) logging if this is turned on
            logAIFile=[obj.filename '_logAI.csv'];
            if obj.aq_settings.acquire_AI
                ch = addAnalogInputChannel(s,device.ID,[0:7], 'Voltage');        
                lh = addlistener(s, 'DataAvailable', @(src, event) obj.logAIData(src, event, logAIFile));
                s.NotifyWhenDataAvailableExceeds = round(s.Rate*1);
                disp(['AI logging enabled. Analog inputs should be connected to ai0 - ai7']);
            end
            
            % Store the session in the FIP_acquisition object
            obj.daq_settings.session=s;
        end
        
        function logAIData(obj, src, ev, logAIFile)
            % LOGAIDATA is responsible for saving analog input (AI) data
            
            dlmwrite(logAIFile,[ev.TimeStamps ev.Data],'delimiter',',','-append');
            
        end
        
        function uncloseable(varargin)
            %UNCLOSABLE makes it impossible to close the window with the
            %stop acquisition button.
            warning('Closing this window will not stop data acquisition. Please press the stop acquisition button instead.')
        end
    end
end
