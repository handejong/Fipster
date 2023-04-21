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
    %       - Photometrics Prime (note: Matlab adapter in Bèta)
    %       - Macvideo facetime camera
    %       - FLIR cameras using Pointgrey adapter (only Flycapture)
    %
    %   FIP_acquisition is part of FIPster. FIPster is made by Johannes de
    %   Jong, j.w.dejong@berkeley.edu
    
    
    properties
        filename        % How the saved datafiles wil be named
        cam_settings    % Struct with camera settings and video source object
        daq_settings    % Struct with DAQ info and session
        aq_settings     % Struct with acquisition settings
        arduino_settings% Struct with arduino_settings
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
            
            % By default we'll use the daq to trigger the setup and camera,
            % but user will be able to change that
            obj.cam_settings.trigger_type = 'daq'; %could be set to arduino
            
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
            
            % Populate Arduino settings
            obj.arduino_settings.port = 'COM3';
            try
                obj.arduino_settings.serial = serial(obj.arduino_settings.port);
            catch
                disp('No Arduino detected, if you intend to use one, change the port under arduino_settings.')
            end
            
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
            obj.plt_settings.type = '-';
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
            
            % if neccesary setup arduino
            if strcmp(obj.cam_settings.trigger_type, 'arduino')
                obj.arduino_setup;
            end
            
            % Where we will store calibration frames
            frames = zeros(res(2), res(1), n_frames);
            
            % Start whatever generates the pulses (DAQ, Arduino or none)
            if strcmp(obj.cam_settings.trigger_type,'daq') %DAQ
                start(obj.cam_settings.vid);
                startBackground(obj.daq_settings.session);
                i=0;
                while i < n_frames
                    i = i + 1;
                    %temp = getdata(obj.cam_settings.vid, 1, 'uint16');
                    frames(:,:,i) = getdata(obj.cam_settings.vid, 1, 'uint16');
                end
                stop(obj.cam_settings.vid);
                stop(obj.daq_settings.session);
                
            elseif strcmp(obj.cam_settings.trigger_type, 'arduino') %Arduino
                start(obj.cam_settings.vid);
                obj.arduino_toggle();
                i = 0;
                while i < n_frames
                    i = i + 1;
                    frames(:,:,i) = getdata(obj.cam_settings.vid, 1, 'uint16');
                end
                stop(obj.cam_settings.vid);
                obj.arduino_toggle(); %Start stop Arduino
                fclose(obj.arduino_settings.serial);
                
            else % No external trigger, just internal
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
            
            % Some cameras will only allow ajusting the ROI by certain
            % intervals. So whe have to check what was actually done.
            ROI = obj.cam_settings.vid.ROIPosition;
            
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
                caxis([50 150])
                colormap gray
            end
            
            % Make a variable were we can put a 'baseline' frame, collected
            % around frame 100, which is used for the life plot and maybe
            % other corrections in the future.
            baseline_frame = uint16(zeros(size(obj.aq_settings.masks(:,:,1))));
            
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
            
            
            % Start the camera
            obj.camera_setup();
            start(obj.cam_settings.vid);
            
            % Start the daq if used as a camera trigger and/or for
            % collection of analogue input data
            if strcmp(obj.cam_settings.trigger_type, 'daq') || obj.aq_settings.acquire_AI
                obj.daq_setup();
                prepare(obj.daq_settings.session); %Reduces start latency
            end
            
            
            % Start the arduino if used as a camera/LED trigger
            if strcmp(obj.cam_settings.trigger_type, 'arduino')
                obj.arduino_setup();
                obj.arduino_toggle();
            end
            
            
            % And start the DAQ at the very last moment. If using both an
            % Arduino and a DAQ they are manually aligned in the Arduino
            % code, there is no trigger. That's not ideal. Should probably
            % also split one of the FIP trigger (LED or camera) into the
            % DAQ to see if your trigger pulses make sense.
            if strcmp(obj.cam_settings.trigger_type, 'daq') || obj.aq_settings.acquire_AI
                startBackground(obj.daq_settings.session)
            end
            
            %%%%%%%%%%%%%%%%%%%%%% Acquisition loop %%%%%%%%%%%%%%%%%%%%%%%
            should_be_time = 0;
            start_time = tic;
            error_counter = 0;
            analysis_time_timer = tic; % only the first time
            while(stop_button.Value == 0)
              
                % Frame number and expected time
                n_frame = n_frame+1;
                should_be_time = should_be_time + 1 / obj.aq_settings.rate;
                
                % Exponentially increase the size of the signal and ref
                % variables.
                if n_frame>length(signal)
                    signal = [signal; zeros(size(signal))];
                    ref=[ref; zeros(size(ref))];
                end
                
                % Get analysis time (loop starts below acquisition)
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
                        
                    case 'arduino'
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
                    disp(['The camera says the last frame is: ' num2str(metadata.FrameNumer) '.'])
                    disp(['The software says the last frame is: ' num2str(n_frame) '.'])    
                end
                
                % This could be the baseline frame
                if n_frame == 100
                    baseline_frame = img - 100;
                end
                
                % Segment the individual fibers out
                for j = 1:obj.aq_settings.fibers
                    data(j) = mean(img(obj.aq_settings.masks(:,:,j)));
                end
                
                % Store the data in the correct channel
                % FOR  TWO-CHANNEL RECORDINGS THE FIRST FRAME IS TAKEN TO
                % BE REF, WHILE THE 2TH FRAME IS THE FRIST SIGNAL FRAME.
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
                    live_feed.CData=img - baseline_frame;
                    
                    % If it's frame 100, callibrate as well
                    if n_frame==100
                        min_value = min(min(live_feed.CData))-20;
                        max_value = max(max(live_feed.CData))+20;
                        figure(live_feed.Parent.Parent)
                        caxis([min_value max_value]);
                    end
                    
                end
                
                % Check if the Arduino (if used) has anything to say
                if strcmp(obj.cam_settings.trigger_type, 'arduino') && obj.arduino_settings.serial.BytesAvailable>0
                    obj.Arduino_serial_reader; % This function will deal with whatever the Arduino has to say
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
            if strcmp(obj.cam_settings.trigger_type,'daq') || obj.aq_settings.acquire_AI
                stop(obj.daq_settings.session)
            end
            
            % If neccesary, stop the arduino and close serial connection
            if strcmp(obj.cam_settings.trigger_type,'arduino')
                obj.arduino_toggle();
                fclose(obj.arduino_settings.serial);
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
            %
            %   Another problem is that the internal clock of the Arduino
            %   is slightly off. For instance, the Arduino MEGA we are
            %   currently using is about 0.2% FASTER than the NI card
            %   (which we assume is perfect). This is a serious delay of
            %   about 7.2sec every hour. To deal with this, if we use both
            %   an Arduino and a DAQ (for logAI scoring) we split the camera
            %   trigger into the DAQ to record exact timetamps of the
            %   camera trigger. We then collect those timestamps.
            
            if strcmp(obj.cam_settings.trigger_type,'arduino')
                % collumn #4 has the camera triggers
                trigger_channel = 4;
                temp_logAI = csvread([obj.filename, '_logAI.csv']);
                indexer = temp_logAI(:,trigger_channel)>1;
                indexer2 = [0; indexer(1:end-1)];
                stamps = temp_logAI(indexer & ~indexer2);
                
                % find the average interval between the stamps
                interval = mean(stamps(2:end)-stamps(1:end-1));
                
                % Update the timeline for signal and/or refference. If
                % there are two channels, the reference channel is the
                % first picture and just the first stamp.
                if obj.aq_settings.channels == 2; signal_start = stamps(2)-0.5*interval; end
                if obj.aq_settings.channels == 1; signal_start = stamps(1)-0.5*interval; end
                ref_start = stamps(1)-0.5 * interval;
                
                % The interval between datapoints is bigger if there are
                % more channels.
                interval = interval*obj.aq_settings.channels;
                
                % Fill out the correct timeline for each signal.
                for i=1:obj.aq_settings.fibers
                    signal(:,2,i)=[signal_start:interval:(length(signal)-1)*interval+signal_start];
                    ref(:,2,i)=[ref_start:interval:(length(ref)-1)*interval+ref_start];
                end
                disp(['Exact timestamps using TTL camera trigger on DAQ collumn: ', num2str(trigger_channel)])
            end
            
            % The start times are 0.5*interval and 1.5 * interval, because the
            % first trigger is at 0, but the first readout is at
            % T = exposure.
            
            % A similar aproach when we used a DAQ to trigger the camera
            % and LEDs
            if strcmp(obj.cam_settings.trigger_type,'daq')
                interval=(1/obj.aq_settings.rate)*obj.aq_settings.channels;
                if obj.aq_settings.channels == 2; signal_start=0.75 * interval; end
                if obj.aq_settings.channels == 1; signal_start = 0.5 * interval; end
                ref_start=0.25 * interval; % Because if there is a ref, there are 2 channels...

                for i=1:obj.aq_settings.fibers
                    signal(:,2,i)=[signal_start:interval:length(signal)*interval];
                    ref(:,2,i)=[ref_start:interval:length(ref)*interval];
                end
                disp('Exact timestamps using  the DAQ')
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
            % Super bad form, but very convenient as well.
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
                case 'pointgrey'
                    supported = true;
                case 'pmimaq_2022b'
                    supported = true;
                case 'gentl'
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
                disp(['Please run calibration before acquisition.'])
            end
   
            % Figure out the correct settings for supported cameras
            switch option.adaptors
                case 'macvideo'
                    % Setup for Apple Facetime cam
                    obj.cam_settings.trigger_type = 'internal';
                    triggerconfig(obj.cam_settings.vid, 'Manual')
                    obj.cam_settings.vid.TriggerRepeat = Inf;
                    obj.cam_settings.vid.FramesPerTrigger=1;
                    obj.cam_settings.vid.ReturnedColorspace = 'grayscale';
                    
                case 'pointgrey'
                    % Setup for a point grey camera
                    
                    % Trigger setup
                    triggerconfig(obj.cam_settings.vid,...
                    'hardware',...
                    'risingEdge',...
                    'externalTriggerMode0-Source0');
                    obj.cam_settings.vid.TriggerRepeat = Inf;
                    obj.cam_settings.vid.FramesPerTrigger = 1;
                    
                    % Exposure settings
                    obj.cam_settings.src.FrameRate = 1;
                    obj.cam_settings.src.FrameRateMode = 'Off';
                    obj.cam_settings.src.ExposureMode = 'Off'; % <- target grey value
                    obj.cam_settings.src.ShutterMode = 'Manual';
                    max_constraint = propinfo(obj.cam_settings.src,'Shutter');
                    max_constraint = max_constraint.ConstraintValue(2); % This is what the camera things is the max <- It's all LIES but I don't know how to change this rignt now
                    obj.cam_settings.src.Shutter = min((1000/obj.aq_settings.rate) - obj.aq_settings.exposure_gap, max_constraint);
                    
                    disp(['Exposure set to: ' num2str(obj.cam_settings.src.Shutter) 'ms']);
                    
                    % Some more settings
                    obj.cam_settings.src.GainMode = 'Manual';
                    obj.cam_settings.src.Gain = 0;
                    try
                        obj.cam_settings.src.SharpnessMode = 'Off';
                        obj.cam_settings.src.GammaMode = 'Off';
                    catch
                        disp('This camera adapter has no Sharpness or Gamma Mode')
                    end

                case 'gentl'
                    % Setup for Gentl (generally Spinakker) cameras
                    
                    % Trigger setup
                    framesPerTrigger = 1;
                    numTriggers = Inf;
                    triggerCondition = "DeviceSpecific";
                    triggerSource = "DeviceSpecific";
                    triggerconfig(obj.cam_settings.vid, "hardware", triggerCondition, triggerSource);
                    obj.cam_settings.vid.FramesPerTrigger = framesPerTrigger;
                    obj.cam_settings.vid.TriggerRepeat = numTriggers - 1;

                    % More trigger settings
                    obj.cam_settings.src.TriggerSelector = "FrameStart";
                    obj.cam_settings.src.TriggerSource = "Line0";
                    obj.cam_settings.src.TriggerActivation = "RisingEdge";
                    obj.cam_settings.src.TriggerMode = "on";

                    % Exposure settings
                    %obj.cam_settings.src.AutoExposureControlLoopDamping = 0;
                    %obj.cam_settings.src.AutoExposureLightingMode = "Backlight";
                    obj.cam_settings.src.AutoExposureTargetGreyValueAuto = "Off";
                    obj.cam_settings.src.DefectCorrectStaticEnable = "False";
                    obj.cam_settings.src.ExposureAuto = "Off";
                    obj.cam_settings.src.GainAuto = "Off";
                    %obj.cam_settings.src.Gamma = 0;
                    obj.cam_settings.src.GammaEnable = "False";
                    obj.cam_settings.src.ExposureMode = "Timed";
                    %obj.cam_settings.src.ExposureTimeAbs = 1000/obj.aq_settings.rate;
                    obj.cam_settings.src.ExposureTime = (1000000/obj.aq_settings.rate) - (1000*obj.aq_settings.exposure_gap); % NOTE micro-seconds

                case 'pmimaq'
                    % Setup for Photometrics
                    
                    % Trigger setup
                    obj.cam_settings.vid.FramesPerTrigger = 1; 
                    obj.cam_settings.vid.TriggerRepeat = Inf;
                    obj.cam_settings.src.TriggerMode='Edge Trigger';
                    
                    % Camera settings
                    obj.cam_settings.src.AutoContrast = 'OFF';
                    obj.cam_settings.src.PP0ENABLED = 'NO';
                    obj.cam_settings.src.PP1ENABLED = 'NO';
                    obj.cam_settings.src.PP2ENABLED = 'NO';
                    obj.cam_settings.src.PP3ENABLED = 'NO';
                    obj.cam_settings.src.PP4ENABLED = 'NO';
                    obj.cam_settings.src.ClearCycles = 1;
                    obj.cam_settings.src.ClearMode = 'Pre-Exposure';
                    
                    % Set the exposure
                    obj.cam_settings.src.Exposure = (1000/obj.aq_settings.rate) - obj.aq_settings.exposure_gap;
                    disp(['Exposure set to: ' num2str(obj.cam_settings.src.Exposure) 'ms']);
          
                    
                    disp('Setup for Prime.')

                case 'pmimaq_2022b'
                    % Setup for Photometrics Prime 2022
                    % THEY CHANGED THE SPELLING OF SOME OF THE PARAMETERS.
                    % WHY?!?!?!
                    
                    % Trigger setup
                    obj.cam_settings.vid.FramesPerTrigger = 1; 
                    obj.cam_settings.vid.TriggerRepeat = Inf;
                    obj.cam_settings.src.TriggerMode='Edge Trigger';
                    
                    % Camera settings
                    triggerconfig(obj.cam_settings.vid, 'hardware')
                    obj.cam_settings.src.AutoContrast = 'OFF';
                    obj.cam_settings.src.PP_0_Enabled = 0;
                    obj.cam_settings.src.PP_1_Enabled = 0;
                    obj.cam_settings.src.PP_2_Enabled = 0;
                    obj.cam_settings.src.PP_3_Enabled = 0;
                    obj.cam_settings.src.PP_4_Enabled = 0;
                    obj.cam_settings.src.ClearCycles = 1;
                    obj.cam_settings.src.ClearMode = 'Pre-Exposure';
                    
                    % Set the exposure
                    obj.cam_settings.src.Exposure = (1000/obj.aq_settings.rate) - obj.aq_settings.exposure_gap;
                    disp(['Exposure set to: ' num2str(obj.cam_settings.src.Exposure) 'ms']);
          
                    
                    disp('Setup for Prime (2022 adapter).')

                otherwise
                    error([option.adaptors ' currently not suported. See cam_setup method.'])
            end
            
            % Update what for which camera the object is currently set
            obj.cam_settings.setup_for=camDeviceN;
        end
        
        function daq_setup(obj)
            %DAQ_SETUP sets upt the daq according to the settings
            
            % Check if there is a daq available
            if strcmp(obj.daq_settings.devices, 'none')
                disp('No DAQ available for AI recording or triggering of the LEDs and Camera.')
                return
            end
            
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
  
            % Setting up the LED triggers (1 or 2 LEDs)
            try %This only works on certain DAQs
                
                % Setting up the camera triggers
                camCh = s.addCounterOutputChannel(device.ID, cam_port, 'PulseGeneration');
                camCh.Frequency = rate;
                camCh.InitialDelay = 0;
                camCh.DutyCycle = 0.1;
                disp(['Camera should be connected to ' camCh.Terminal]);
                
                % Setting up LED triggers
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
            catch
                warning('This DAQ can not be used for camera or LED triggering, but AI data is recorded.')
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
            obj.daq_settings.session = s;
        end
        
        function arduino_setup(obj)
            % Sets up the arduino
            
            % NOTE. In oder for this to work, the accompanying Arduino code
            % needs to be loaded on an Arduino. See the notes inside that
            % Arduino sketch for more info on the serial communication
            % protocol.
            
            % Open serial connection
            s = obj.arduino_settings.serial;
            fopen(s);
            pause(1);
            
            % Check if proper program running...
            
%             % Set tech delay
%             t_delay = obj.aq_settings.exposure_gap;
%             fwrite(s,3); fwrite(s,t_delay); fgets(s) %3 is the tech delay noun
            
            % Set framerate
            rate=obj.aq_settings.rate;
            fwrite(s,2); fwrite(s,rate); fgets(s) %2 is the framerate noun
            
        end
        
        function logAIData(obj, src, ev, logAIFile)
            % LOGAIDATA is responsible for saving analog input (AI) data
            
            dlmwrite(logAIFile,[ev.TimeStamps ev.Data],'delimiter',',','-append');
            
        end
        
        function arduino_toggle(obj)
            % Will turn on and off arduino sampling
            
            s = obj.arduino_settings.serial;
            
            % Check if serial is properly working
            if s.bytesavailable>0
                error('Serial connection error')
            end
            
            % 9 is Acquisition, 1 is empty verb or stop if running
            fwrite(s,9); fwrite(s,1); fgets(s)

        end
        
        function Arduino_serial_reader(obj)
            % This function is called when there are bytes available on the
            % serial connection between the Arduino and Matlab. It
            % processes whatever the Arduino has to say.
            
            disp('The Arduino says: ');
            disp(fgets(obj.arduino_settings.serial)); % Just prints the data
            
        end
        
        function uncloseable(varargin)
            %UNCLOSABLE makes it impossible to close the window with the
            %stop acquisition button.
            warning('Closing this window will not stop data acquisition. Please press the stop acquisition button instead.')
        end
    end
end

