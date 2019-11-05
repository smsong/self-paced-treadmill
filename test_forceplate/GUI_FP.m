function varargout = GUI_FP(varargin)
% GUI_FP MATLAB code for GUI_FP.fig
%      GUI_FP, by itself, creates a new GUI_FP or raises the existing
%      singleton*.
%
%      H = GUI_FP returns the handle to a new GUI_FP or the handle to
%      the existing singleton*.
%
%      GUI_FP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FP.M with the given input arguments.
%
%      GUI_FP('Property','Value',...) creates a new GUI_FP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_FP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_FP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_FP

% Last Modified by GUIDE v2.5 04-Nov-2019 18:27:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_FP_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_FP_OutputFcn, ...
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


% --- Executes just before GUI_FP is made visible.
function GUI_FP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_FP (see VARARGIN)
global flag_GRFCalibrated saved_grf_offsets

% Choose default command line output for GUI_FP
handles.output = hObject;

% Speedgoat setup
startopt = questdlg('Start options for Self-Paced Treadmill','Start Options', ...
    'Build','Load','Start','Start');
switch startopt
    case 'Build'
        slbuild('ForcePlate');
        start(slrt);
    case 'Load'
        load(slrt,'ForcePlate');
        start(slrt);
    otherwise
        start(slrt);
end

handles.tg = SimulinkRealTime.target;

setappdata(hObject, 'flag_plotFP_R', 1);
setappdata(hObject, 'flag_plotFP_L', 1);

handles.dt_timerPlot = .05;
handles.dt_timerGRFCalibrate = .05;
handles.timer_plots = timer('TimerFcn',{@timerPlots, hObject}, ...
    'Period', handles.dt_timerPlot , ...
    'ExecutionMode','fixedRate','BusyMode','queue');
handles.timer_GRFCalibrate = timer('TimerFcn',{@timerGRFCalibrate, hObject}, ...
    'Period', handles.dt_timerGRFCalibrate , ...
    'ExecutionMode','fixedRate','BusyMode','queue');

if flag_GRFCalibrated
    setappdata(hObject, 'grf_offsets', saved_grf_offsets);
	set(handles.btn_tm_calibrate, 'BackgroundColor', 'green');
else
    setappdata(hObject, 'grf_offsets', zeros(1,12));
end
setappdata(hObject, 'n_grfs_for_cal', 20);
setappdata(hObject, 'i_grfs_for_cal', 0);

scale_grf = [-500, -500, 1000, -800, -400, -400]; % V -> N
setappdata(hObject, 'scale_grf', scale_grf);

% Update handles structure
guidata(hObject, handles);

% start timers
start(handles.timer_plots);

% UIWAIT makes GUI_FP wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_FP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% =============== %
% TIMER FUNCTIONS %
% =============== %
% timerPlots(...)
% timerGRFCalibrate(...)

function timerPlots(~, ~, hObject)
if getappdata(hObject, 'flag_plotFP_R')
    updateFPRPlot(hObject);
end
if getappdata(hObject, 'flag_plotFP_L')
    updateFPLPlot(hObject);
end


function handles = timerGRFCalibrate(~, ~, hObject)
global flag_GRFCalibrated saved_grf_offsets
handles = guidata(hObject);

n_grf = getappdata(hObject, 'n_grfs_for_cal');
i_grf = getappdata(hObject, 'i_grfs_for_cal');

fxl_raw = getsignal(handles.tg,'sensor_TML/s1');
fyl_raw = getsignal(handles.tg,'sensor_TML/s2');
fzl_raw = getsignal(handles.tg,'sensor_TML/s3');
mxl_raw = getsignal(handles.tg,'sensor_TML/s4');
myl_raw = getsignal(handles.tg,'sensor_TML/s5');
mzl_raw = getsignal(handles.tg,'sensor_TML/s6');
fxr_raw = getsignal(handles.tg,'sensor_TMR/s1');
fyr_raw = getsignal(handles.tg,'sensor_TMR/s2');
fzr_raw = getsignal(handles.tg,'sensor_TMR/s3');
mxr_raw = getsignal(handles.tg,'sensor_TMR/s4');
myr_raw = getsignal(handles.tg,'sensor_TMR/s5');
mzr_raw = getsignal(handles.tg,'sensor_TMR/s6');
    
grfs = [fxl_raw fyl_raw fzl_raw mxl_raw myl_raw mzl_raw ...
    fxr_raw fyr_raw fzr_raw mxr_raw myr_raw mzr_raw];

grf_offsets = getappdata(hObject, 'grf_offsets');
grf_offsets = (grf_offsets*i_grf + grfs)/(i_grf + 1);
setappdata(hObject, 'grf_offsets', grf_offsets);

i_grf = i_grf + 1;

if i_grf >= n_grf
    i_grf = 0;
    set(handles.btn_tm_calibrate, 'BackgroundColor', 'green');
    stop(handles.timer_GRFCalibrate);
    flag_GRFCalibrated = 1;
    saved_grf_offsets = grf_offsets;
end

setappdata(hObject, 'i_grfs_for_cal', i_grf);
guidata(hObject, handles);



% ============= %
% GUI FUNCTIONS %
% ============= %
% updateFPRPlot(...)

function updateFPRPlot(hObject)
handles = guidata(hObject);
LEG = 1; % -1: L; 1: R

fx_raw = getsignal(handles.tg,'sensor_TMR/s1');
fy_raw = getsignal(handles.tg,'sensor_TMR/s2');
fz_raw = getsignal(handles.tg,'sensor_TMR/s3');
mx_raw = getsignal(handles.tg,'sensor_TMR/s4');
my_raw = getsignal(handles.tg,'sensor_TMR/s5');
mz_raw = getsignal(handles.tg,'sensor_TMR/s6');

scale_grf = getappdata(hObject, 'scale_grf');
grf_offsets = getappdata(hObject, 'grf_offsets');

grfs(1) = scale_grf(1)*(fx_raw - grf_offsets(6+1));
grfs(2) = scale_grf(2)*(fy_raw - grf_offsets(6+2));
grfs(3) = scale_grf(3)*(fz_raw - grf_offsets(6+3));
grfs(4) = scale_grf(4)*(mx_raw - grf_offsets(6+4));
grfs(5) = scale_grf(5)*(my_raw - grf_offsets(6+5));
grfs(6) = scale_grf(6)*(mz_raw - grf_offsets(6+6));

h_axes = handles.axes_FP_R;
data = get(handles.axes_FP_R, 'UserData');
data = updateFPPlot(h_axes, LEG, data, grfs);

set(handles.axes_FP_R, 'UserData', data);
guidata(hObject, handles);

function updateFPLPlot(hObject)
handles = guidata(hObject);
LEG = -1; % -1: L; 1: R

fx_raw = getsignal(handles.tg,'sensor_TML/s1');
fy_raw = getsignal(handles.tg,'sensor_TML/s2');
fz_raw = getsignal(handles.tg,'sensor_TML/s3');
mx_raw = getsignal(handles.tg,'sensor_TML/s4');
my_raw = getsignal(handles.tg,'sensor_TML/s5');
mz_raw = getsignal(handles.tg,'sensor_TMR/s6');

scale_grf = getappdata(hObject, 'scale_grf');
grf_offsets = getappdata(hObject, 'grf_offsets');

grfs(1) = scale_grf(1)*(fx_raw - grf_offsets(1));
grfs(2) = scale_grf(2)*(fy_raw - grf_offsets(2));
grfs(3) = scale_grf(3)*(fz_raw - grf_offsets(3));
grfs(4) = scale_grf(4)*(mx_raw - grf_offsets(4));
grfs(5) = scale_grf(5)*(my_raw - grf_offsets(5));
grfs(6) = scale_grf(6)*(mz_raw - grf_offsets(6));

h_axes = handles.axes_FP_L;
data = get(handles.axes_FP_L, 'UserData');
data = updateFPPlot(h_axes, LEG, data, grfs);

set(handles.axes_FP_L, 'UserData', data);
guidata(hObject, handles);



function data = updateFPPlot(h_axes, LEG, data, grfs)

h_fp = 0.015; % forceplate height [m]
plate_xy = [0.5588; % plate width = 22 inch
            1.8034]; % plate length = 71 inch
cop_off = [ .5*plate_xy(1); % plate width = 22 inch
            -.5*plate_xy(2)]; % plate length = 71 inch

if ~isfield(data, 'flag_init')
    hold(h_axes, 'on');
    h_data_cop = plot(h_axes, ...
                0, 0, 'k.', ...
                'Markersize', 20);
    xlim(h_axes, .5*plate_xy(1)*[-1 1]);
    ylim(h_axes, .5*plate_xy(2)*[-1 1]);
    
    data.flag_init = 1;
    data.h_data_cop = h_data_cop;
end

fx = grfs(1);
fy = grfs(2);
fz = grfs(3);
mx = grfs(4);
my = grfs(5);

cop(1) = -(h_fp*fx + my)/fz + LEG*cop_off(1);
cop(2) = -(h_fp*fy + mx)/fz + cop_off(2);

if fz > 10 % [N]
    data.h_data_cop.YData = cop(2);
    data.h_data_cop.XData = cop(1);
    data.h_data_cop.Color = [1 0 0];
else
    data.h_data_cop.Color = [1 1 1];
end



% ========== %
% GUI INPUTS %
% ========== %

function btn_abort_Callback(hObject, eventdata, handles)
stop(handles.timer_plots);
delete(handles.output);


function btn_tm_calibrate_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
set(handles.btn_tm_calibrate, 'BackgroundColor', 'yellow');
start(handles.timer_GRFCalibrate);
