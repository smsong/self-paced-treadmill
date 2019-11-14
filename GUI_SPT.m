function varargout = GUI_SPT(varargin)
% GUI_SPT MATLAB code for GUI_SPT.fig
%      GUI_SPT, by itself, creates a new GUI_SPT or raises the existing
%      singleton*.
%
%      H = GUI_SPT returns the handle to a new GUI_SPT or the handle to
%      the existing singleton*.
%
%      GUI_SPT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SPT.M with the given input arguments.
%
%      GUI_SPT('Property','Value',...) creates a new GUI_SPT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_SPT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_SPT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_SPT

% Last Modified by GUIDE v2.5 04-Nov-2019 12:52:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_SPT_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_SPT_OutputFcn, ...
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


% --- Executes just before GUI_SPT is made visible.
function GUI_SPT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_SPT (see VARARGIN)
global flag_GRFCalibrated flag_subMass saved_grf_offsets saved_mass ...
    flag_ctrl_par ctrl_par

% Choose default command line output for GUI_SPT
handles.output = hObject;

% Speedgoat setup
startopt = questdlg('Start options for Self-Paced Treadmill','Start Options', ...
    'Build','Load','Start','Start');
switch startopt
    case 'Build'
        slbuild('Treadmill_SPT');
        start(slrt);
    case 'Load'
        load(slrt,'Treadmill_SPT');
        start(slrt);
    otherwise
        start(slrt);
end

format compact
handles.tg = SimulinkRealTime.target;

setappdata(hObject, 'n_reset', 0);

setappdata(hObject, 'flag_SPT', 0); % SPT: self-paced treadmill

setappdata(hObject, 'datetime_start', clock);

handles.dt_timerPlot = .05;
handles.dt_timerTreadmillCtrl = .01;
handles.dt_timerGRFCalibrate = .05;
handles.timer_plots = timer('TimerFcn',{@timerPlots, hObject}, ...
    'Period', handles.dt_timerPlot , ...
    'ExecutionMode','fixedRate','BusyMode','queue');
handles.timer_treadmillCtrl = timer('TimerFcn',{@timerTreadmillCtrl, hObject}, ...
    'Period', handles.dt_timerTreadmillCtrl , ...
    'ExecutionMode','fixedRate','BusyMode','queue');
handles.timer_GRFCalibrate = timer('TimerFcn',{@timerGRFCalibrate, hObject}, ...
    'Period', handles.dt_timerGRFCalibrate , ...
    'ExecutionMode','fixedRate','BusyMode','queue');
handles.timer_subMass = timer('TimerFcn',{@timerSubMass, hObject}, ...
    'Period', handles.dt_timerGRFCalibrate , ...
    'ExecutionMode','fixedRate','BusyMode','queue');

if flag_GRFCalibrated
    setappdata(hObject, 'grf_offsets', saved_grf_offsets);
	set(handles.btn_tm_calibrate, 'BackgroundColor', 'green');
    setparam(handles.tg, 'Treadmill_Ctrl/grf_offsets', 'Value', saved_grf_offsets);
else
    setappdata(hObject, 'grf_offsets', zeros(1,12));
end
setappdata(hObject, 'n_grfs_for_cal', 20);
setappdata(hObject, 'i_grfs_for_cal', 0);
setappdata(hObject, 'n_grfs_for_mass', 20);
setappdata(hObject, 'i_grf_for_mass', 0);

if flag_ctrl_par
    setparam(handles.tg, 'Treadmill_Ctrl/p_off', 'Value', ctrl_par(1));    
    setparam(handles.tg, 'Treadmill_Ctrl/G_p', 'Value', ctrl_par(2));
    setparam(handles.tg, 'Treadmill_Ctrl/G_v', 'Value', ctrl_par(3));
    setparam(handles.tg, 'Treadmill_Ctrl/r_v_tm_tgt', 'Value', [ctrl_par(4), ctrl_par(5)]);
    setparam(handles.tg, 'Treadmill_Ctrl/del_t_tgt', 'Value', ctrl_par(6));

    set(handles.edit_ctrl_p_off, 'String', ctrl_par(1));
    set(handles.edit_ctrl_G_p, 'String', ctrl_par(2));
    set(handles.edit_ctrl_G_v, 'String', ctrl_par(3));
    set(handles.edit_ctrl_min_v, 'String', ctrl_par(4));
    set(handles.edit_ctrl_max_v, 'String', ctrl_par(5));
    set(handles.edit_ctrl_del_t, 'String', ctrl_par(6));
else    
    p_off0 = 0;
    G_p0 = 0.1;
    G_v0 = 0.25;
    min_v0 = 0;
    max_v0 = 3.0;
    del_t0 = 0.5;
    
    ctrl_par(1) = p_off0;
    ctrl_par(2) = G_p0;
    ctrl_par(3) = G_v0;
    ctrl_par(4) = min_v0;
    ctrl_par(5) = max_v0;
    ctrl_par(6) = del_t0;
end

if flag_subMass
    mass = saved_mass;
    set(handles.btn_tm_measure_mass, 'BackgroundColor', 'green');
    setparam(handles.tg, 'GUI/mass', 'Value', mass);
else
    mass = 80; % [kg]; only for initialization
end
setappdata(hObject, 'mass', mass);
g = 9.81;
setappdata(hObject, 'g', g); % [m/s^2]; gravity

scale_grf = [-500, -500, 1000, -800, -400, -400]; % V -> N
setappdata(hObject, 'scale_grf', scale_grf); 

setappdata(hObject, 'scale_fy_plot', scale_grf(2)/(mass*g));
setappdata(hObject, 'scale_fz_plot', scale_grf(3)/(mass*g));

mass = getappdata(hObject, 'mass');
scale_grf = getappdata(hObject, 'scale_grf');
scale_fy = scale_grf(2)/mass;
scale_fz = scale_grf(3)/mass;
scale_mx = scale_grf(4)/mass;
scale_fy_fz_mx = [scale_fy, scale_fz, scale_mx];
setappdata(hObject, 'scale_fy_fz_mx', scale_fy_fz_mx);
setparam(handles.tg, 'Treadmill_Ctrl/scale_fy_fz_mx', 'Value', scale_fy_fz_mx);

% setappdata(hObject, 'm2mm', 1000); % m -> mm
% setappdata(hObject, 'mm2m', .001); % m -> mm

setappdata(hObject, 'time0_ctrl', 0);
setappdata(hObject, 'v_tm', 0); % [m/s] treadmill speed
setappdata(hObject, 'v_tm_tgt', 0); % [m/s] target treadmill speed
setappdata(hObject, 'a_tm0', .5); % [m/s^2] treadmill accelerate (used to change speed)

% file scope (to save data)
handles.fscope = getscope(handles.tg, 101);
set(handles.fscope, 'NumSamples', 60*60*1000);
handles.filename = 'TM_test.dat';
set(handles.fscope, 'filename', handles.filename);

% Update handles structure
guidata(hObject, handles);

% start timers
start(handles.timer_treadmillCtrl);
start(handles.timer_plots);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_SPT_OutputFcn(hObject, eventdata, handles) 
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
% timerTreadmillCtrl(...)
% timerGRFCalibrate(...)
% timerSubMass(...)

function timerPlots(~, ~, hObject)

handles = guidata(hObject);
if getappdata(hObject, 'flag_SPT')
    datetime_start = getappdata(hObject, 'datetime_start_spt');
    str_main = 'Self Paced Mode';
    updateBTNDisplay(handles.btn_tm_spt, datetime_start, str_main);
        
    datetime = clock;
    time_now = datetime(4)*60*60 + datetime(5)*60 + datetime(6);
        
    datetime_start_10 = getappdata(hObject, 'datetime_start_spt_del_10');
    time_start_del_10 = datetime_start_10(4)*60*60 + datetime_start_10(5)*60 + datetime_start_10(6);
    del_t_10 = time_now - time_start_del_10;
    if del_t_10 >= 10
        p_tm0 = getappdata(handles.output, 'x_10sec');
        p_tm = getsignal(handles.tg,'Treadmill_Ctrl/p_tm/s1');        
        
        fprintf('average v (10 sec): %.3f [m/s]\n', (p_tm - p_tm0)/del_t_10);
        
        setappdata(hObject, 'x_10sec', p_tm);
        setappdata(hObject, 'datetime_start_spt_del_10', datetime);
    end
    
    datetime_start_60 = getappdata(hObject, 'datetime_start_spt_del_60');
    time_start_del_60 = datetime_start_60(4)*60*60 + datetime_start_60(5)*60 + datetime_start_60(6);
    del_t_60 = time_now - time_start_del_60;
    if del_t_60 >= 60
        p_tm0 = getappdata(handles.output, 'x_60sec');
        p_tm = getsignal(handles.tg,'Treadmill_Ctrl/p_tm/s1');        
        
        fprintf('average v (60 sec): %.3f [m/s]\n', (p_tm - p_tm0)/del_t_60);
        
        setappdata(hObject, 'x_60sec', p_tm);
        setappdata(hObject, 'datetime_start_spt_del_60', datetime);
    end
end

drawnow;


function timerTreadmillCtrl(~, ~, hObject)
handles = guidata(hObject);

if getappdata(hObject, 'flag_SPT')
    v_tm_tgt = getsignal(handles.tg,'Treadmill_Ctrl/v_tm_tgt_out/s1');
    v_tm_tgt0 = getappdata(hObject, 'v_tm_tgt');

    if  v_tm_tgt ~= v_tm_tgt0
        a_tm_tgt = getsignal(handles.tg,'Treadmill_Ctrl/a_tm_tgt_out/s1');

        setTreadmillSpeed(handles, v_tm_tgt, a_tm_tgt);

        setappdata(hObject, 'v_tm_tgt', v_tm_tgt);
    end
end

function handles = timerGRFCalibrate(~, ~, hObject)
global flag_GRFCalibrated saved_grf_offsets
handles = guidata(hObject);

n_grf = getappdata(hObject, 'n_grfs_for_cal');
i_grf = getappdata(hObject, 'i_grfs_for_cal');

fxl_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TML/s1');
fyl_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TML/s2');
fzl_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TML/s3');
mxl_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TML/s4');
myl_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TML/s5');
mzl_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TML/s6');
fxr_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TMR/s1');
fyr_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TMR/s2');
fzr_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TMR/s3');
mxr_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TMR/s4');
myr_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TMR/s5');
mzr_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TMR/s6');
    
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
    setparam(handles.tg, 'Treadmill_Ctrl/grf_offsets', 'Value', saved_grf_offsets);
end

setappdata(hObject, 'i_grfs_for_cal', i_grf);
guidata(hObject, handles);


function handles = timerSubMass(~, ~, hObject)
global flag_subMass saved_mass
handles = guidata(hObject);

n_grf = getappdata(hObject, 'n_grfs_for_mass');
i_grf = getappdata(hObject, 'i_grf_for_mass');

g = getappdata(hObject, 'g');

scale_grf = getappdata(hObject, 'scale_grf');
scale_fz__ = scale_grf(3);

fzl_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TML/s3');
fzr_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TMR/s3');

grf_offsets = getappdata(hObject, 'grf_offsets');
mass1 = scale_fz__/g*(fzl_raw - grf_offsets(3) ...
        + fzr_raw - grf_offsets(6+3));

mass = getappdata(hObject, 'mass');
mass = (mass*i_grf + mass1)/(i_grf + 1);
setappdata(hObject, 'mass', mass);

i_grf = i_grf + 1;

if i_grf >= n_grf
    i_grf = 0;
    disp(num2str(mass));
    set(handles.btn_tm_measure_mass, 'BackgroundColor', 'green');
    stop(handles.timer_subMass);    
    
    scale_grf = getappdata(hObject, 'scale_grf');
    scale_fy = scale_grf(2)/mass;
    scale_fz = scale_grf(3)/mass;
    scale_mx = scale_grf(4)/mass;
    
    scale_fy_fz_mx = [scale_fy, scale_fz, scale_mx];    
    setappdata(hObject, 'scale_fy_fz_mx', scale_fy_fz_mx);
    setparam(handles.tg, 'Treadmill_Ctrl/scale_fy_fz_mx', 'Value', scale_fy_fz_mx);
    
    setappdata(hObject, 'scale_fy_plot', scale_grf(2)/(mass*g));
    setappdata(hObject, 'scale_fz_plot', scale_grf(3)/(mass*g));
    
    flag_subMass = 1;
    saved_mass = mass;
    setparam(handles.tg, 'GUI/mass', 'Value', mass);
end

setappdata(hObject, 'i_grf_for_mass', i_grf);
guidata(hObject, handles);


% ============= %
% GUI FUNCTIONS %
% ============= %
% updateBTNDisplay(...)


function updateBTNDisplay(handle_btn, datetime_start, str_main)
datetime = clock;
time_start_spt = datetime_start(4)*60*60 + datetime_start(5)*60 + datetime_start(6);
time_now = datetime(4)*60*60 + datetime(5)*60 + datetime(6);
time_spt = round(time_now - time_start_spt);
time_spt_hms(1) = floor(time_spt/60/60); % hour
time_spt_m = time_spt-time_spt_hms(1)*60*60;
time_spt_hms(2) = floor(time_spt_m/60); % minutes
time_spt_s = time_spt_m - time_spt_hms(2)*60;
time_spt_hms(3) = time_spt_s; % hour

s_hms = {'', '', ''};
for i = 1:3
    s_hms{i} = num2str(time_spt_hms(i));
    if length(s_hms{i}) == 1
        s_hms{i} = ['0' s_hms{i}];
    end
end
display_gui = ['<html><div style="text-align:center">' str_main '<br>' ...
    s_hms{1} ':' s_hms{2} ':' s_hms{3} '</div></html>'];
set(handle_btn, 'String', display_gui)


% ============== %
% USER FUNCTIONS %
% ============== %
% setTreadmillSpeed(...)


function setTreadmillSpeed(handles, v_tgt, a_tgt)
global ctrl_par
min_v = ctrl_par(4);
max_v = ctrl_par(5);
v_tgt = max(min_v, min(max_v, v_tgt));

% this function does not work in high frequency
m2mm = 1000;

defInc = 0; % WARNING: DO NOT CHANGE THIS UNLESS YOU KNOW WHAT YOU ARE DOING

speedR = v_tgt*m2mm;
speedL = v_tgt*m2mm;
accR = a_tgt*m2mm;
accL = a_tgt*m2mm;
incline = defInc;

h_tcpip = openTreadmillComm();
[payload] = getPayload(speedR, speedL, accR, accL, incline);
sendTreadmillPacket(payload, h_tcpip);
closeTreadmillComm(h_tcpip);

% send new treadmill speed to speedgoat
setparam(handles.tg, 'Treadmill_Ctrl/v_tm_tgt_in', 'Value', v_tgt);
setparam(handles.tg, 'Treadmill_Ctrl/a_tm_tgt_in', 'Value', a_tgt);

set(handles.edit_tm_v, 'String', round(v_tgt, 2));


function startSelfPacedMode(handles)
set(handles.btn_tm_spt, 'BackgroundColor', 'green');
datetime = clock;
setappdata(handles.output, 'datetime_start_spt', datetime);
setappdata(handles.output, 'datetime_start_spt_del_10', datetime);
setappdata(handles.output, 'datetime_start_spt_del_60', datetime);

p_tm = getsignal(handles.tg,'Treadmill_Ctrl/p_tm/s1');
setappdata(handles.output, 'x_10sec', p_tm);
setappdata(handles.output, 'x_60sec', p_tm);

flag_SPT = 1;    
set(handles.btn_tm_set_v,'Enable','off')
setappdata(handles.output, 'flag_SPT', flag_SPT); 
setparam(handles.tg, 'GUI/flag_SPT', 'Value', flag_SPT);


function stopSelfPacedMode(handles)
set(handles.btn_tm_spt, 'BackgroundColor', [0.94, 0.94, 0.94]);
set(handles.btn_tm_set_v,'Enable','on')

flag_SPT = 0;    
setappdata(handles.output, 'flag_SPT', flag_SPT); 
setparam(handles.tg, 'GUI/flag_SPT', 'Value', flag_SPT);



% =========== %
% GUI BUTTONS %
% =========== %

% --- Executes on button press in btn_tm_set_v.
function btn_tm_set_v_Callback(hObject, eventdata, handles)
str_v_tgt = get(handles.edit_tm_tgt_v, 'String');
%check if the input is a number. if so, send command to treadmill
if( ~isnan(str2double(str_v_tgt)))
    v_tgt = str2double(str_v_tgt);
    a_tm0 = getappdata(handles.output, 'a_tm0');
    setTreadmillSpeed(handles, v_tgt, a_tm0);
    setappdata(handles.output, 'v_tm_tgt', v_tgt);
else
    disp('target velocity should be a number!!!');
end

guidata(hObject, handles); % !!!

% --- Executes on button press in btn_tm_stop.
function btn_tm_stop_Callback(hObject, eventdata, handles)
flag_SPT = getappdata(handles.output, 'flag_SPT');
if flag_SPT == 1 
    stopSelfPacedMode(handles)
end
a_tm0 = getappdata(handles.output, 'a_tm0');
setTreadmillSpeed(handles, 0, a_tm0);
setappdata(handles.output, 'v_tm_tgt', 0);

% --- Executes on button press in btn_tm_calibrate.
function btn_tm_calibrate_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
set(handles.btn_tm_calibrate, 'BackgroundColor', 'yellow');
start(handles.timer_GRFCalibrate);

% --- Executes on button press in btn_KF_reset.
function btn_KF_reset_Callback(hObject, eventdata, handles)
n_reset = getappdata(handles.output, 'n_reset');
n_reset = n_reset + 1;
setappdata(handles.output, 'n_reset', n_reset);

setparam(handles.tg, 'Treadmill_Ctrl/trackCOM_reset_trigger', 'Value', n_reset);

% --- Executes on button press in btn_tm_measure_mass.
function btn_tm_measure_mass_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
set(handles.btn_tm_measure_mass, 'BackgroundColor', 'yellow');
start(handles.timer_subMass);

% --- Executes on button press in btn_tm_spt.
function btn_tm_spt_Callback(hObject, eventdata, handles)
flag_SPT = getappdata(handles.output, 'flag_SPT');
if flag_SPT == 0
    startSelfPacedMode(handles);
elseif flag_SPT == 1
    stopSelfPacedMode(handles);
end

% --- Executes on button press in btn_abort.
function btn_abort_Callback(hObject, eventdata, handles)
try
    % stop all timers
    if strcmp(handles.timer_treadmillCtrl.Running, 'on')
        stop(handles.timer_treadmillCtrl);
    end
    if strcmp(handles.timer_plots.Running, 'on')
        stop(handles.timer_plots);
    end
    if strcmp(handles.timer_GRFCalibrate.Running, 'on')
        stop(handles.timer_GRFCalibrate);
    end
    if strcmp(handles.timer_subMass.Running, 'on')
        stop(handles.timer_subMass);
    end

    flag_SPT = getappdata(handles.output, 'flag_SPT');
    if flag_SPT == 1 
        stopSelfPacedMode(handles)
    end

    % stop and close treadmill
    a_tm0 = getappdata(handles.output, 'a_tm0');
    setTreadmillSpeed(handles, 0, a_tm0);
    setappdata(handles.output, 'v_tm_tgt', 0);
catch exception
    disp('gui not closed properly');
    disp(exception)
end
delete(handles.output);


function edit_tm_tgt_v_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tm_tgt_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tm_tgt_v as text
%        str2double(get(hObject,'String')) returns contents of edit_tm_tgt_v as a double


% --- Executes during object creation, after setting all properties.
function edit_tm_tgt_v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tm_tgt_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_tm_v_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tm_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tm_v as text
%        str2double(get(hObject,'String')) returns contents of edit_tm_v as a double


% --- Executes during object creation, after setting all properties.
function edit_tm_v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tm_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% h = findobj('Tag', 'btn_abort');
btn_abort_Callback(handles.btn_abort, eventdata, handles)

% Hint: delete(hObject) closes the figure
delete(hObject);


function edit_ctrl_G_p_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ctrl_G_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ctrl_G_p as text
%        str2double(get(hObject,'String')) returns contents of edit_ctrl_G_p as a double


% --- Executes during object creation, after setting all properties.
function edit_ctrl_G_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ctrl_G_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_ctrl_set.
function btn_ctrl_set_Callback(hObject, eventdata, handles)
p_off_new = str2double( get(handles.edit_ctrl_p_off, 'String') );
G_p_new = str2double( get(handles.edit_ctrl_G_p, 'String') );
G_v_new = str2double( get(handles.edit_ctrl_G_v, 'String') );
min_v_new = str2double( get(handles.edit_ctrl_min_v, 'String') );
max_v_new = str2double( get(handles.edit_ctrl_max_v, 'String') );
del_t_new = str2double( get(handles.edit_ctrl_del_t, 'String') );

setparam(handles.tg, 'Treadmill_Ctrl/p_off', 'Value', p_off_new);
setparam(handles.tg, 'Treadmill_Ctrl/G_p', 'Value', G_p_new);
setparam(handles.tg, 'Treadmill_Ctrl/G_v', 'Value', G_v_new);
setparam(handles.tg, 'Treadmill_Ctrl/r_v_tm_tgt', 'Value', [min_v_new, max_v_new]);
setparam(handles.tg, 'Treadmill_Ctrl/del_t_tgt', 'Value', del_t_new);

global flag_ctrl_par ctrl_par
ctrl_par(1) = p_off_new;
ctrl_par(2) = G_p_new;
ctrl_par(3) = G_v_new;
ctrl_par(4) = min_v_new;
ctrl_par(5) = max_v_new;
ctrl_par(6) = del_t_new;
flag_ctrl_par = 1;


function edit_ctrl_G_v_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ctrl_G_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ctrl_G_v as text
%        str2double(get(hObject,'String')) returns contents of edit_ctrl_G_v as a double


% --- Executes during object creation, after setting all properties.
function edit_ctrl_G_v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ctrl_G_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_ctrl_min_v_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ctrl_min_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ctrl_min_v as text
%        str2double(get(hObject,'String')) returns contents of edit_ctrl_min_v as a double


% --- Executes during object creation, after setting all properties.
function edit_ctrl_min_v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ctrl_min_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_ctrl_max_v_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ctrl_max_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ctrl_max_v as text
%        str2double(get(hObject,'String')) returns contents of edit_ctrl_max_v as a double


% --- Executes during object creation, after setting all properties.
function edit_ctrl_max_v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ctrl_max_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_ctrl_del_t_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ctrl_del_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ctrl_del_t as text
%        str2double(get(hObject,'String')) returns contents of edit_ctrl_del_t as a double


% --- Executes during object creation, after setting all properties.
function edit_ctrl_del_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ctrl_del_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_ctrl_reset.
function btn_ctrl_reset_Callback(hObject, eventdata, handles)
p_off0 = 0;
G_p0 = 0.1;
G_v0 = 0.25;
min_v0 = 0;
max_v0 = 3.0;
del_t0 = 0.5;

setparam(handles.tg, 'Treadmill_Ctrl/p_off', 'Value', p_off0);
setparam(handles.tg, 'Treadmill_Ctrl/G_p', 'Value', G_p0);
setparam(handles.tg, 'Treadmill_Ctrl/G_v', 'Value', G_v0);
setparam(handles.tg, 'Treadmill_Ctrl/r_v_tm_tgt', 'Value', [min_v0, max_v0]);
setparam(handles.tg, 'Treadmill_Ctrl/del_t_tgt', 'Value', del_t0);

set(handles.edit_ctrl_p_off, 'String', p_off0);
set(handles.edit_ctrl_G_p, 'String', G_p0);
set(handles.edit_ctrl_G_v, 'String', G_v0);
set(handles.edit_ctrl_min_v, 'String', min_v0);
set(handles.edit_ctrl_max_v, 'String', max_v0);
set(handles.edit_ctrl_del_t, 'String', del_t0);

global flag_ctrl_par ctrl_par
ctrl_par(1) = p_off0;
ctrl_par(2) = G_p0;
ctrl_par(3) = G_v0;
ctrl_par(4) = min_v0;
ctrl_par(5) = max_v0;
ctrl_par(6) = del_t0;
flag_ctrl_par = 1;


function edit_ctrl_p_off_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ctrl_p_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ctrl_p_off as text
%        str2double(get(hObject,'String')) returns contents of edit_ctrl_p_off as a double


% --- Executes during object creation, after setting all properties.
function edit_ctrl_p_off_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ctrl_p_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
