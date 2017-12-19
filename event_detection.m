function varargout = event_detection(varargin)
% EVENT_DETECTION MATLAB code for event_detection.fig
%      EVENT_DETECTION, by itself, creates a new EVENT_DETECTION or raises the existing
%      singleton*.
%
%      H = EVENT_DETECTION returns the handle to a new EVENT_DETECTION or the handle to
%      the existing singleton*.
%
%      EVENT_DETECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVENT_DETECTION.M with the given input arguments.
%
%      EVENT_DETECTION('Property','Value',...) creates a new EVENT_DETECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before event_detection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to event_detection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help event_detection

% Last Modified by GUIDE v2.5 23-Jan-2017 14:51:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @event_detection_OpeningFcn, ...
                   'gui_OutputFcn',  @event_detection_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before event_detection is made visible.
function event_detection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to event_detection (see VARARGIN)

% Choose default command line output for event_detection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Storing the handle to the sweepset object and select it
handles.paired_sweepset=varargin{1};
figure(handles.paired_sweepset.handles.figure)

% Drawing all elements, but not making them visible
x_values=xlim;
y_values=ylim;
%   -Baseline
%   -Event markers


% Add listener
handles.listener(1)=addlistener(handles.paired_sweepset,'state_change',@(scr, ev) update_sweep(scr, ev, handles));
handles.listener(2)=addlistener(handles.paired_sweepset,'selection_change',@(scr, ev) update_everything(scr, ev, handles));
handles.listener(3)=addlistener(handles.paired_sweepset,'baseline_change',@(scr, ev) update_everything(scr, ev, handles));

% Update handles structure
guidata(hObject, handles);

% Filling variable in the GUI
set(handles.filename,'String',handles.paired_sweepset.filename)

% Find events
event_list=find_events(handles);
setappdata(hObject,'event_list',event_list);

% Setting a callback for then this GUI is closed and other callbacks
set(hObject,'CloseRequestFcn',{@close_req, handles})


% --- Outputs from this function are returned to the command line.
function varargout = event_detection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function Threshold_slider_Callback(hObject, eventdata, handles)
% hObject    handle to Threshold_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function Threshold_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Threshold_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% just to make sure that not only this GUI, but also it's associated
% doodles are closed.

function close_req(src,ev,handles)
delete(handles.display_handles.event_markers)
delete(handles.display_handles.baseline)

 
delete(src)
delete(handles.listener)
 
 
function update_everything(scr,ev,handles)

find_peak(handles);

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Other %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function find_events(handles)
paired_sweepset=handles.paired_sweepset;

   
SF=paired_sweepset.sampling_frequency;





set(handles.display_handles.peak_up,'Position',[peak.location_up,peak.up,0],'String',display_text_up);
set(handles.display_handles.peak_down,'Position',[peak.location_down,peak.down,0],'String',display_text_down);

% updating displayed vallues
set(handles.maximum,'String',num2str(peak.up))
set(handles.minimum,'String',num2str(peak.down))
