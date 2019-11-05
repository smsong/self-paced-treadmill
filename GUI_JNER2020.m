function varargout = GUI_JNER2020(varargin)
% GUI_JNER2020 MATLAB code for GUI_JNER2020.fig
%      GUI_JNER2020, by itself, creates a new GUI_JNER2020 or raises the existing
%      singleton*.
%
%      H = GUI_JNER2020 returns the handle to a new GUI_JNER2020 or the handle to
%      the existing singleton*.
%
%      GUI_JNER2020('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_JNER2020.M with the given input arguments.
%
%      GUI_JNER2020('Property','Value',...) creates a new GUI_JNER2020 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_JNER2020_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_JNER2020_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_JNER2020

% Last Modified by GUIDE v2.5 04-Nov-2019 14:40:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_JNER2020_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_JNER2020_OutputFcn, ...
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


% --- Executes just before GUI_JNER2020 is made visible.
function GUI_JNER2020_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_JNER2020 (see VARARGIN)
global flag_GRFCalibrated flag_subMass saved_grf_offsets saved_mass ...
    flag_ctrl_par ctrl_par

% Choose default command line output for GUI_JNER2020
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

setappdata(hObject, 'v0', 0.8);

setappdata(hObject, 'n_reset', 0);

setappdata(hObject, 'flag_plotCom', 0);
setappdata(hObject, 'com', 0);
setappdata(hObject, 'dcom', 0);
setappdata(hObject, 'y_step', 0);
setappdata(hObject, 'leg_step', 'L'); % 'L' or 'R'

setappdata(hObject, 'flag_SPT', 0); % SPT: self-paced treadmill
setappdata(hObject, 'flag_sweep', 0); % treadmill speed sweep mode
setappdata(hObject, 'flag_tgt_distance', 0); % display target distance
setappdata(hObject, 'flag_tgt_time', 0); % display target distance
setappdata(hObject, 'flag_1min_started', 0);
setappdata(hObject, 'flag_100m_started', 0);
setappdata(hObject, 'flag_mean_p', 0); % track mean values
setappdata(hObject, 'i_mark_TMSP2', 0);
setappdata(hObject, 't_mark_TMSP2', 0);
setappdata(hObject, 'i_mark_TMSP150', 0);
setappdata(hObject, 't_mark_TMSP150', 0);

setappdata(hObject, 'mean_p_n', 0);
setappdata(hObject, 'mean_p_p', 0);

setappdata(hObject, 'del_t_tgt', .5);

setappdata(hObject, 'datetime_start', clock);

setappdata(handles.output, 'tgtDist_track', 150);
setappdata(handles.output, 'tgtDist_laps', 1);

setappdata(hObject, 'n_step', 0);

setappdata(hObject, 'x_10sec', 0);
setappdata(hObject, 'x_60sec', 0);

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
    set(handles.axes_COM_KF, 'Color', 'white');
    setappdata(hObject, 'flag_plotCom', 1);
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

setappdata(hObject, 'time0_ctrl', 0);
setappdata(hObject, 'v_tm', 0); % [m/s] treadmill speed
setappdata(hObject, 'v_tm_tgt', 0); % [m/s] target treadmill speed
setappdata(hObject, 'a_tm0', .5); % [m/s^2] treadmill accelerate (used to change speed)

setappdata(hObject, 'flag_saveData', 0); %flag for saving data
setappdata(hObject, 'save_index', 1);

% file scope (to save data)
handles.fscope = getscope(handles.tg, 101);
set(handles.fscope, 'NumSamples', 60*60*1000);
handles.filename = 'TM_test.dat';
set(handles.fscope, 'filename', handles.filename);
handles.data_saved = {};

% Update handles structure
guidata(hObject, handles);

% start timers
start(handles.timer_treadmillCtrl);
start(handles.timer_plots);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_JNER2020_OutputFcn(hObject, eventdata, handles) 
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
updateLGRFPlots(hObject);

if getappdata(hObject, 'flag_plotCom')
    updateSubjectPlots(hObject);
end

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

if getappdata(hObject, 'flag_tgt_time')
    datetime_start = getappdata(hObject, 'datetime_start_tgtTime');
    str_main = 'SP 2 min';
    updateBTNDisplay(handles.btn_SP2, datetime_start, str_main);
    
    time_start = datetime_start(4)*60*60 + datetime_start(5)*60 + datetime_start(6);
    datetime = clock;
    time_now = datetime(4)*60*60 + datetime(5)*60 + datetime(6);
    t_tm_track = time_now - time_start;
    
    tgtTime_i_lap = getappdata(handles.output, 'tgtTime_i_lap');    
    p_tm = getsignal(handles.tg,'Treadmill_Ctrl/p_tm/s1');
    p_tm0 = getappdata(handles.output, 'p_mark_TMSP2');
    
    flag_1min_started = getappdata(hObject, 'flag_1min_started');
    if ~flag_1min_started && t_tm_track > (tgtTime_i_lap)*120 + 60
        setappdata(hObject, 'p_1min', p_tm);
        setappdata(hObject, 'flag_1min_started', 1);
    end
    v_t_mark = [6 14 20 26 34 40 46 54 60 66 74 80 86 94 100 106 114 120];
    i_mark_TMSP2 = getappdata(hObject, 'i_mark_TMSP2');

    datetime = clock;
    time_now = datetime(4)*60*60 + datetime(5)*60 + datetime(6);
    t_mark_TMSP2 = getappdata(hObject, 't_mark_TMSP2');
    if time_now - t_mark_TMSP2 > (tgtTime_i_lap)*120 + v_t_mark(i_mark_TMSP2+1)
        i_mark_TMSP2 = i_mark_TMSP2 + 1;
        fprintf('#%d: %.3f [m]\n', ...
            i_mark_TMSP2, p_tm - p_tm0);
        if i_mark_TMSP2 >= length(v_t_mark)
            i_mark_TMSP2 = i_mark_TMSP2 - length(v_t_mark);
        end
        setappdata(hObject, 'i_mark_TMSP2', i_mark_TMSP2);
        setappdata(hObject, 'p_mark_TMSP2', p_tm);    
    end
    if time_now - t_mark_TMSP2 > (tgtTime_i_lap+1)*120
        tgtTime_i_lap = tgtTime_i_lap + 1;

        p_start_tgtTime_lab = getappdata(hObject, 'p_start_tgtTime_lab');
        del_p = p_tm-p_start_tgtTime_lab;
        v_time = del_p/120;

        p_start_1min = getappdata(hObject, 'p_1min');
        del_p_1min = p_tm - p_start_1min;
        v_time_1min = del_p_1min/60;

        datetime = clock;
    
        fprintf('lab %d: %.1f [m], %.3f [m/s], [last 1 min] %.3f [m/s]\n', ...
            tgtTime_i_lap, del_p, v_time, v_time_1min);
        setappdata(hObject, 'tgtTime_i_lap', tgtTime_i_lap);
        setappdata(hObject, 'datetime_start_tgtTime_lab', datetime);
        setappdata(hObject, 'flag_1min_started', 0);
    end
end

if getappdata(hObject, 'flag_tgt_distance')
    datetime_start = getappdata(hObject, 'datetime_start_tgtDist');
    str_main = 'SP 150 m';
    updateBTNDisplay(handles.btn_SP150, datetime_start, str_main);
    
    if ishandle(handles.fig_tgt_distance)
        updateTgtDistancePlot(hObject);
    end
end

if getappdata(hObject, 'flag_sweep')
    datetime_start = getappdata(hObject, 'datetime_start_sweep');
    str_main = 'MSS';
    updateBTNDisplay(handles.btn_tm_sweep_v, datetime_start, str_main);
end

if getappdata(hObject, 'flag_mean_p')
    datetime_start = getappdata(hObject, 'datetime_start_mean');
    str_main = 'Mean p_off';
    updateBTNDisplay(handles.btn_mean_p, datetime_start, str_main);
end

drawnow;


function timerTreadmillCtrl(~, ~, hObject)
handles = guidata(hObject);

if getappdata(hObject, 'flag_SPT')
    % flag_step = getsignal(handles.tg,'Treadmill_Ctrl/flag_step/s1');
    % can set a trigger from speedgoat?
    v_tm_tgt = getsignal(handles.tg,'Treadmill_Ctrl/v_tm_tgt_out/s1');
    v_tm_tgt0 = getappdata(hObject, 'v_tm_tgt');

    if  v_tm_tgt ~= v_tm_tgt0
        a_tm_tgt = getsignal(handles.tg,'Treadmill_Ctrl/a_tm_tgt_out/s1');

        setTreadmillSpeed(handles, v_tm_tgt, a_tm_tgt);

        setappdata(hObject, 'v_tm_tgt', v_tm_tgt);
    end
end

if getappdata(hObject, 'flag_sweep')
    del_t_sweep = str2double(get(handles.edit_tm_del_t, 'String'));    
    datetime_start = getappdata(hObject, 'datetime_start_sweep_del');
    datetime = clock;
    time_start_del = datetime_start(4)*60*60 + datetime_start(5)*60 + datetime_start(6);
    time_now = datetime(4)*60*60 + datetime(5)*60 + datetime(6);
    del_t = time_now - time_start_del;
    if del_t >= del_t_sweep
        v_tm = str2double(get(handles.edit_tm_v, 'String'));
        del_v = str2double(get(handles.edit_tm_del_v, 'String'));
        setTreadmillSpeed(handles, v_tm + del_v, .5);
        
        setappdata(hObject, 'datetime_start_sweep_del', datetime);
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
    
    % start com tracking
    if strcmp(handles.timer_treadmillCtrl.Running, 'off')
        set(handles.axes_COM_KF, 'Color', 'white');
        start(handles.timer_treadmillCtrl);
        setappdata(hObject, 'flag_plotCom', 1);
    end
end

setappdata(hObject, 'i_grf_for_mass', i_grf);
guidata(hObject, handles);


% ============= %
% GUI FUNCTIONS %
% ============= %
% updateLGRFPlots(...)
% updateSubjectPlots(...)
% updateBTNDisplay(...)

function updateLGRFPlots(hObject)
handles = guidata(hObject);
data = get(handles.axes_GRF_L, 'UserData');

t_min = 0;
t_max = 2.0;
data_len = round((t_max-t_min)/handles.dt_timerPlot);

if ~isfield(data, 'flag_init')
    data_grfy = [linspace(t_min, t_max, data_len);
                                nan(1, data_len)];
    data_grfz = [linspace(t_min, t_max, data_len);
                                nan(1, data_len)];
    hold(handles.axes_GRF_L, 'on');
    h_data_grfy = plot(handles.axes_GRF_L, ...
        data_grfy(1,:), data_grfy(2,:), ...
        'r-', 'LineWidth', 1);
    h_data_grfz = plot(handles.axes_GRF_L, ...
        data_grfz(1,:), data_grfz(2,:), ...
        'b-', 'LineWidth', 1);
    ylim(handles.axes_GRF_L, [-.5 2]);
    xlim(handles.axes_GRF_L, [t_min t_max]);
    
    data.flag_init = 1;
    data.h_data_grfy = h_data_grfy;
    data.h_data_grfz = h_data_grfz;
    data.time0 = 0;
    data.n_step = 0;
    data.i = 1;

    mass = getappdata(hObject, 'mass');
    g = getappdata(hObject, 'g');
end

time = getsignal(handles.tg,'Treadmill_Data_Export/sensor_time/s1');

if data.i > length(data.h_data_grfy.XData)
    data.time0 = time;
    data.i = 1;
    data.h_data_grfy.YData = nan(1, data_len);
    data.h_data_grfz.YData = nan(1, data_len);
end

fyl_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TML/s2');
fzl_raw = getsignal(handles.tg,'Treadmill_Data_Export/sensor_TML/s3');

grf_offsets = getappdata(hObject, 'grf_offsets');
scale_fy_plot = getappdata(hObject, 'scale_fy_plot');
scale_fz_plot = getappdata(hObject, 'scale_fz_plot');
fyl_ = scale_fy_plot*(fyl_raw - grf_offsets(2));
fzl_ = scale_fz_plot*(fzl_raw - grf_offsets(3));

data.h_data_grfy.XData(data.i) = time - data.time0;
data.h_data_grfy.YData(data.i) = fyl_;
data.h_data_grfz.XData(data.i) = time - data.time0;
data.h_data_grfz.YData(data.i) = fzl_;

data.i = data.i + 1;

set(handles.axes_GRF_L, 'UserData', data);
guidata(hObject, handles);


function updateSubjectPlots(hObject)
handles = guidata(hObject);
data = get(handles.axes_COM_KF, 'UserData');

if ~isfield(data, 'flag_init')
    hold(handles.axes_COM_KF, 'on');
    h_data_com = plot(handles.axes_COM_KF, ...
                0, 0, 'k.', ...
                'Markersize', 20);
    h_data_dcom = quiver(handles.axes_COM_KF, ...
                0, 0, 0, 0, 'k-', ...
                'LineWidth', 1);
    h_data_p_step_L = plot(handles.axes_COM_KF, ...
                -.1, 0, 'r.', ...
                'Markersize', 20);
    h_data_p_step_R = plot(handles.axes_COM_KF, ...
                .1, 0, 'b.', ...
                'Markersize', 20);
    plot(handles.axes_COM_KF, [-.15 .15], [0 0], 'k:');
    ylim(handles.axes_COM_KF, [-1 1]);
    xlim(handles.axes_COM_KF, [-.15 .15]);
    
    data.flag_init = 1;
    data.h_data_com = h_data_com;
    data.h_data_dcom = h_data_dcom;
    data.h_data_p_step_L = h_data_p_step_L;
    data.h_data_p_step_R = h_data_p_step_R;
end

com = getsignal(handles.tg,'Treadmill_Ctrl/com/s1');
dcom = getsignal(handles.tg,'Treadmill_Ctrl/dcom/s1');
y_step_1 = getsignal(handles.tg,'Treadmill_Ctrl/y_step_1/s1');
y_step_0 = getsignal(handles.tg,'Treadmill_Ctrl/y_step_0/s1');
leg_step = getsignal(handles.tg,'Treadmill_Ctrl/leg_step/s1');
stance_L = getsignal(handles.tg,'Treadmill_Ctrl/stance_L/s1');
stance_R = getsignal(handles.tg,'Treadmill_Ctrl/stance_R/s1');

data.h_data_com.YData = com;
data.h_data_dcom.YData = com;
data.h_data_dcom.VData = dcom;
if leg_step % 0: L; 1: R
    y_step_R = y_step_1;
    y_step_L = y_step_0;
else
    y_step_L = y_step_1;
    y_step_R = y_step_0;
end

if stance_L
    data.h_data_p_step_L.YData = y_step_L;
    data.h_data_p_step_L.Color = [1 0 0];
else
    data.h_data_p_step_L.Color = [1 1 1];
end
if stance_R
    data.h_data_p_step_R.YData = y_step_R;
    data.h_data_p_step_R.Color = [0 0 1];
else
    data.h_data_p_step_R.Color = [1 1 1];
end


if getappdata(hObject, 'flag_mean_p')
    mean_p_n = getappdata(hObject, 'mean_p_n');
    mean_p_p = getappdata(hObject, 'mean_p_p');
    
    mean_p_p = (mean_p_p*mean_p_n + com)/(mean_p_n + 1);
    mean_p_n = mean_p_n + 1;
    
    setappdata(hObject, 'mean_p_n', mean_p_n);
    setappdata(hObject, 'mean_p_p', mean_p_p);
end

set(handles.axes_COM_KF, 'UserData', data);
guidata(hObject, handles);

function updateTgtDistancePlot(hObject)
handles = guidata(hObject);
p_tm = getsignal(handles.tg,'Treadmill_Ctrl/p_tm/s1');
p_tm0 = getappdata(handles.output, 'tgtDist_p_tm_0');
tgtDist_track = getappdata(handles.output, 'tgtDist_track');

p_tm_track = p_tm - p_tm0;
p_tm_track_i = mod(p_tm_track, tgtDist_track);
track_l = 100*tgtDist_track/400;
track_r = track_l/pi;

if p_tm_track_i < tgtDist_track/8
    handles.h_tgtDist_pos.XData = p_tm_track_i;
    handles.h_tgtDist_pos.YData = track_r;
elseif p_tm_track_i < tgtDist_track*3/8
    handles.h_tgtDist_pos.XData = track_l/2 + track_r*sin((p_tm_track_i-tgtDist_track/8)*4*pi/tgtDist_track);
    handles.h_tgtDist_pos.YData = track_r*cos((p_tm_track_i-tgtDist_track/8)*4*pi/tgtDist_track);
elseif p_tm_track_i < tgtDist_track*5/8
    handles.h_tgtDist_pos.XData = track_l/2 - p_tm_track_i + tgtDist_track*3/8;
    handles.h_tgtDist_pos.YData = -track_r;
elseif p_tm_track_i < tgtDist_track*7/8
    handles.h_tgtDist_pos.XData = -track_l/2 + track_r*sin((p_tm_track_i-tgtDist_track*5/8)*4*pi/tgtDist_track+pi);
    handles.h_tgtDist_pos.YData = track_r*cos((p_tm_track_i-tgtDist_track*5/8)*4*pi/tgtDist_track+pi);
else
    handles.h_tgtDist_pos.XData = -track_l/2 + p_tm_track_i - tgtDist_track*7/8;
    handles.h_tgtDist_pos.YData = track_r;
end

tgtDist_i_lap = getappdata(handles.output, 'tgtDist_i_lap');
% print velocity during last 100 m
flag_100m_started = getappdata(hObject, 'flag_100m_started');
if ~flag_100m_started && p_tm_track > (tgtDist_i_lap)*tgtDist_track + 50
    datetime = clock;
    setappdata(hObject, 'datetime_50m', datetime);
    setappdata(hObject, 'flag_100m_started', 1);
end

v_d_mark = [7.5 17.5 25 32.5 42.5 50 57.5 67.5 75 82.5 92.5 100 107.5 117.5 125 132.5 142.5 150];
i_mark_TMSP150 = getappdata(hObject, 'i_mark_TMSP150');
if p_tm_track > (tgtDist_i_lap)*tgtDist_track + v_d_mark(i_mark_TMSP150+1)
    i_mark_TMSP150 = i_mark_TMSP150 + 1;
    t_mark_TMSP150 = getappdata(hObject, 't_mark_TMSP150');    
    datetime = clock;
    time_now = datetime(4)*60*60 + datetime(5)*60 + datetime(6);
    fprintf('#%d: %.3f [sec]\n', ...
        i_mark_TMSP150, time_now - t_mark_TMSP150);
    if i_mark_TMSP150 >= length(v_d_mark)
        i_mark_TMSP150 = i_mark_TMSP150 - length(v_d_mark);
    end
    setappdata(hObject, 'i_mark_TMSP150', i_mark_TMSP150);
    setappdata(hObject, 't_mark_TMSP150', time_now);
end

if p_tm_track > (tgtDist_i_lap+1)*tgtDist_track
    tgtDist_i_lap = tgtDist_i_lap + 1;
    
    datetime_start = getappdata(hObject, 'datetime_start_tgtDist_lab');
    datetime = clock;
    time_start_distance = datetime_start(4)*60*60 + datetime_start(5)*60 + datetime_start(6);
    time_now = datetime(4)*60*60 + datetime(5)*60 + datetime(6);
    del_t = time_now-time_start_distance;
    v_track = tgtDist_track/del_t;
    
    datetime_start_100m = getappdata(hObject, 'datetime_50m');
    time_start_100m = datetime_start_100m(4)*60*60 + datetime_start_100m(5)*60 + datetime_start_100m(6);
    del_t_100m = time_now-time_start_100m;
    v_track_100m = 100/del_t_100m;
    
    fprintf('lab %d: %.1f [sec], %.3f [m/s], [last 100 m] %.3f [m/s]\n', ...
        tgtDist_i_lap, del_t, v_track, v_track_100m);
    setparam(handles.tg, 'GUI/VT_i_lap', 'Value', tgtDist_i_lap);
    setparam(handles.tg, 'GUI/VT_del_t', 'Value', del_t);
    setparam(handles.tg, 'GUI/VT_v', 'Value', v_track);
    tgtDist_laps = getappdata(handles.output, 'tgtDist_laps');
    handles.h_tgtDist_txt.String = [num2str(tgtDist_i_lap) '/' num2str(tgtDist_laps)];
    setappdata(hObject, 'tgtDist_i_lap', tgtDist_i_lap);
    setappdata(hObject, 'datetime_start_tgtDist_lab', datetime);
    setappdata(hObject, 'flag_100m_started', 0);
end
guidata(hObject, handles);


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
% trackTreadmillSpeed(...)
% setTreadmillSpeed(...)

% only call in timerTreadmillCtrl (uses 'time0_ctrl')
function v_tm = trackTreadmillSpeed(hObject, time)
v_tm = getappdata(hObject, 'v_tm');
v_tm_tgt = getappdata(hObject, 'v_tm_tgt');
a_tm0 = getappdata(hObject, 'a_tm0');
time0 = getappdata(hObject, 'time0_ctrl');

del_t = time - time0;
del_v = v_tm_tgt - v_tm;
if del_v > 0
    del_v = min(del_v, a_tm0*del_t);
else
    del_v = max(del_v, -a_tm0*del_t);
end
v_tm = v_tm + del_v;
setappdata(hObject, 'v_tm', v_tm);


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


function startManualSpeedSelection(handles)
set(handles.btn_tm_sweep_v, 'BackgroundColor', 'green');
datetime = clock;
setappdata(handles.output, 'datetime_start_sweep', datetime);
setappdata(handles.output, 'datetime_start_sweep_del', datetime);

set(handles.btn_tm_set_v,'Enable','off')
    
flag_sweep = 1;
setappdata(handles.output, 'flag_sweep', flag_sweep);
setparam(handles.tg, 'GUI/flag_sweep', 'Value', flag_sweep);


function stopManualSpeedSelection(handles)
set(handles.btn_tm_sweep_v, 'BackgroundColor', [0.94, 0.94, 0.94]);
v_tm = getsignal(handles.tg,'Treadmill_Ctrl/v_tm/s1');
fprintf('final speed: %.3f [m/s]\n', v_tm);
setparam(handles.tg, 'GUI/SS_final_speed', 'Value', v_tm);

% set next starting velosity and del_v
del_v = str2double(get(handles.edit_tm_del_v, 'String'));
set(handles.edit_tm_tgt_v, 'String', num2str(v_tm + sign(del_v)*0.1));
set(handles.edit_tm_del_v, 'String', num2str(-del_v));

set(handles.btn_tm_set_v,'Enable','on')

flag_sweep = 0;
setappdata(handles.output, 'flag_sweep', flag_sweep);
setparam(handles.tg, 'GUI/flag_sweep', 'Value', flag_sweep);


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

flag_tgt_time = getappdata(handles.output, 'flag_tgt_time');
if flag_tgt_time == 1
    flag_tgt_time = 0;
    set(handles.btn_SP2, 'BackgroundColor', [0.94, 0.94, 0.94]);

    setappdata(handles.output, 'flag_tgt_time', flag_tgt_time);
    setparam(handles.tg, 'GUI/flag_TMSP2', 'Value', flag_tgt_time);
    guidata(hObject, handles);
    
    flag_SPT = getappdata(handles.output, 'flag_SPT');
    if flag_SPT == 1
        stopSelfPacedMode(handles);
    end
end

flag_tgt_distance = getappdata(handles.output, 'flag_tgt_distance');
if flag_tgt_distance == 1
    flag_tgt_distance = 0;
    set(handles.btn_SP150, 'BackgroundColor', [0.94, 0.94, 0.94]);
    
    if ishandle(handles.fig_tgt_distance)
        close(handles.fig_tgt_distance)
    end
    
    setappdata(handles.output, 'flag_tgt_distance', flag_tgt_distance);
    setparam(handles.tg, 'GUI/flag_track', 'Value', flag_tgt_distance);
    guidata(hObject, handles);
    
    flag_SPT = getappdata(handles.output, 'flag_SPT');
    if flag_SPT == 1 
        stopSelfPacedMode(handles)
    end
end
    
flag_SPT = getappdata(handles.output, 'flag_SPT');
if flag_SPT == 1 
    stopSelfPacedMode(handles)
end

flag_sweep = getappdata(handles.output, 'flag_sweep');
if flag_sweep == 1
    stopManualSpeedSelection(handles);
end
a_tm0 = getappdata(handles.output, 'a_tm0');
setTreadmillSpeed(handles, 0, a_tm0);
setappdata(handles.output, 'v_tm_tgt', 0);

v0 = getappdata(handles.output, 'v0');
set(handles.edit_tm_tgt_v, 'String', num2str(v0));

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

if strcmp(handles.timer_treadmillCtrl.Running, 'off')
    start(handles.timer_treadmillCtrl);
end

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
guidata(hObject, handles);

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
    
    if isfield(handles, 'fig_tgt_distance') && ishandle(handles.fig_tgt_distance)
        close(handles.fig_tgt_distance)
    end

    flag_SPT = getappdata(handles.output, 'flag_SPT');
    if flag_SPT == 1
        stopSelfPacedMode(handles)
    end
    
    flag_sweep = getappdata(handles.output, 'flag_sweep');
    if flag_sweep == 1
        stopManualSpeedSelection(handles);
    end

    % stop and close treadmill
    a_tm0 = getappdata(handles.output, 'a_tm0');
    setTreadmillSpeed(handles, 0, a_tm0);
    setappdata(handles.output, 'v_tm_tgt', 0);

    
    % stop saving data (if not stopped)
    flag_saveData = getappdata(handles.output, 'flag_saveData');
    if flag_saveData == 1
        stop(handles.fscope);
    end
    
    % export data
    data_saved = unique(handles.data_saved);
    data_path0 = cd;
    if ~isempty(data_saved)
        datetime = clock;
        subdir = ['data_exp\' num2str(datetime(2)) '_' num2str(datetime(3))];
        data_path = [data_path0 '\' subdir];
        try
            cd(data_path);
        catch
            mkdir(subdir);
            cd(data_path);
        end
    end
    for i=1:length(data_saved)%handles.scopes)
        try
            file_server = ['c:\' data_saved{i}];
            SimulinkRealTime.copyFileToHost(handles.tg, file_server);
            disp([data_saved{i} ' copied']);
        catch exception
            disp([data_saved{i} ' not copied']);
            disp(exception)
        end
    end
    cd(data_path0);
catch exception
    disp('gui not closed properly: check data export');
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


% --- Executes on button press in btn_save_data.
function btn_save_data_Callback(hObject, eventdata, handles)
% hObject    handle to btn_save_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flag_saveData = getappdata(handles.output, 'flag_saveData');

if flag_saveData == 0
    % start saving data
    filename_save = get(handles.edit_subject_code, 'String');
    handles.filename = [filename_save '.dat'];
    set(handles.fscope, 'filename', handles.filename);
    start(handles.fscope);
    handles.data_saved = [handles.data_saved {handles.filename}];
    
    % update GUI
    datetime = clock;
    handles.fscope_start = datetime(4:6);
    handles.fscope_disp = {handles.filename; ...
        ['Start: ' num2str(datetime(4)) ':' num2str(datetime(5)) '.' num2str(datetime(6))]};
    set(handles.text_dataFileInfo, 'String', handles.fscope_disp);
    
    setappdata(handles.output, 'flag_saveData', 1);
    set(handles.btn_save_data, 'BackgroundColor', 'green');
    
elseif flag_saveData == 1
    % stop saving data
    stop(handles.fscope);
    
    % update GUI
    datetime = clock;
    handles.fscope_disp = [handles.fscope_disp; ['Stop: ' num2str(datetime(4)) ':' num2str(datetime(5)) '.' num2str(datetime(6))]];
    set(handles.text_dataFileInfo, 'String', handles.fscope_disp);
    setappdata(handles.output, 'flag_saveData', 0);  
    set(handles.btn_save_data, 'BackgroundColor', 'yellow');
end

guidata(hObject, handles);


function btn_flag_save_Callback(hObject, eventdata, handles)
% hObject    handle to btn_flag_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of btn_flag_save as text
%        str2double(get(hObject,'String')) returns contents of btn_flag_save as a double


function edit_subject_code_Callback(hObject, eventdata, handles)
% hObject    handle to edit_subject_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_subject_code as text
%        str2double(get(hObject,'String')) returns contents of edit_subject_code as a double


% --- Executes during object creation, after setting all properties.
function edit_subject_code_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_subject_code (see GCBO)
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


% --- Executes during object deletion, before destroying properties.
function text_dataFileInfo_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to text_dataFileInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% h = findobj('Tag', 'btn_abort');
btn_abort_Callback(handles.btn_abort, eventdata, handles)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in btn_SP2.
function btn_SP2_Callback(hObject, eventdata, handles)
flag_tgt_time = getappdata(handles.output, 'flag_tgt_time');
if flag_tgt_time == 0
    flag_tgt_time = 1;
    set(handles.btn_SP2, 'BackgroundColor', 'green');
    datetime = clock;
    time_now = datetime(4)*60*60 + datetime(5)*60 + datetime(6);
    setappdata(handles.output, 't_mark_TMSP2', time_now);
    setappdata(handles.output, 'datetime_start_tgtTime', datetime);
    
    p_tm = getsignal(handles.tg,'Treadmill_Ctrl/p_tm/s1');
    setappdata(handles.output, 'p_mark_TMSP2', p_tm);
    setappdata(handles.output, 'p_start_tgtTime_lab', p_tm);
    setappdata(handles.output, 'tgtTime_i_lap', 0);

    flag_SPT = getappdata(handles.output, 'flag_SPT');
    if flag_SPT == 0
        startSelfPacedMode(handles);
    end

elseif flag_tgt_time == 1
    flag_tgt_time = 0;
    set(handles.btn_SP2, 'BackgroundColor', [0.94, 0.94, 0.94]);
    
    flag_SPT = getappdata(handles.output, 'flag_SPT');
    if flag_SPT == 1
        stopSelfPacedMode(handles);
    end
end

setappdata(handles.output, 'flag_tgt_time', flag_tgt_time);
setparam(handles.tg, 'GUI/flag_TMSP2', 'Value', flag_tgt_time);
guidata(hObject, handles);


% --- Executes on button press in btn_SP150.
function btn_SP150_Callback(hObject, eventdata, handles)
flag_tgt_distance = getappdata(handles.output, 'flag_tgt_distance');
if flag_tgt_distance == 0    
    % set distance display
    str_track = get(handles.edit_track, 'String');
    %check if the input is a number. if so, send command to treadmill
    if(~isnan(str2double(str_track)))
        tgtDist_track = str2double(str_track);
    else
        disp('track [m] should be a number!!!');
        return;
    end
    str_labs = get(handles.edit_labs, 'String');
    %check if the input is a number. if so, send command to treadmill
    if(~isnan(str2double(str_labs)))
        tgtDist_laps = str2double(str_labs);
    else
        disp('labs should be a number!!!');
        return;
    end
    
    % plot track
    handles.fig_tgt_distance = figure('Position', [100 100 1000 600], 'NumberTitle', 'off');
    set(handles.fig_tgt_distance, 'MenuBar', 'none');
    set(handles.fig_tgt_distance, 'ToolBar', 'none');    
    set(handles.fig_tgt_distance, 'Color', 'w');
    
    track_l = 100*tgtDist_track/400;
    track_r = track_l/pi;
    th1 = 0:pi/50:pi;
    th2 = pi:pi/50:2*pi;
    x_track = [linspace(-track_l/2, track_l/2, 51) track_l/2+track_r*sin(th1) linspace(track_l/2, -track_l/2, 51) -track_l/2+track_r*sin(th2)];
    y_track = [track_r*ones(1,51) track_r*cos(th1) -track_r*ones(1,51) track_r*cos(th2)];
    plot(x_track, y_track, 'm-', 'LineWidth', 10);
    axis off;
    axis equal;
    hold on;
    plot([0 0], track_r*[.8 1.2], 'k-', 'LineWidth', 3);
    p_tm = getsignal(handles.tg,'Treadmill_Ctrl/p_tm/s1');
    handles.h_tgtDist_pos = plot(0, track_r, 'k.', 'MarkerSize', 100);
    handles.h_tgtDist_txt = text(0, 0, ['0/' num2str(tgtDist_laps)], ...
        'FontSize', 60, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
    setappdata(handles.output, 'tgtDist_track', tgtDist_track);
    setappdata(handles.output, 'tgtDist_laps', tgtDist_laps);
    setappdata(handles.output, 'tgtDist_p_tm_0', p_tm);
    setappdata(handles.output, 'tgtDist_i_lap', 0);
        
    flag_tgt_distance = 1;
    set(handles.btn_SP150, 'BackgroundColor', 'green');
    datetime = clock;
    time_now = datetime(4)*60*60 + datetime(5)*60 + datetime(6);
    setappdata(handles.output, 't_mark_TMSP150', time_now);
    setappdata(handles.output, 'datetime_start_tgtDist', datetime);
    setappdata(handles.output, 'datetime_start_tgtDist_lab', datetime);

    flag_SPT = getappdata(handles.output, 'flag_SPT');
    if flag_SPT == 0
        startSelfPacedMode(handles)
    end
elseif flag_tgt_distance == 1
    flag_tgt_distance = 0;
    set(handles.btn_SP150, 'BackgroundColor', [0.94, 0.94, 0.94]);
    
    if ishandle(handles.fig_tgt_distance)
        close(handles.fig_tgt_distance)
    end
    
    flag_SPT = getappdata(handles.output, 'flag_SPT');
    if flag_SPT == 1 
        stopSelfPacedMode(handles)
    end
end
setappdata(handles.output, 'flag_tgt_distance', flag_tgt_distance);
setparam(handles.tg, 'GUI/flag_track', 'Value', flag_tgt_distance);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function axes_COM_KF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_COM_KF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_COM_KF

set(hObject, 'XTick',[], 'YTick',[]);


% --- Executes during object creation, after setting all properties.
function axes_GRF_L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_GRF_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_GRF_L

set(hObject, 'XTick',[], 'YTick',[]);


function edit_track_Callback(hObject, eventdata, handles)
% hObject    handle to edit_track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_track as text
%        str2double(get(hObject,'String')) returns contents of edit_track as a double


% --- Executes during object creation, after setting all properties.
function edit_track_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_labs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_labs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_labs as text
%        str2double(get(hObject,'String')) returns contents of edit_labs as a double


% --- Executes during object creation, after setting all properties.
function edit_labs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_labs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_tm_del_v_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tm_del_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tm_del_v as text
%        str2double(get(hObject,'String')) returns contents of edit_tm_del_v as a double


% --- Executes during object creation, after setting all properties.
function edit_tm_del_v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tm_del_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tm_del_t_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tm_del_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tm_del_t as text
%        str2double(get(hObject,'String')) returns contents of edit_tm_del_t as a double


% --- Executes during object creation, after setting all properties.
function edit_tm_del_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tm_del_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_tm_sweep_v.
function btn_tm_sweep_v_Callback(hObject, eventdata, handles)
flag_sweep = getappdata(handles.output, 'flag_sweep');
if flag_sweep == 0
    startManualSpeedSelection(handles);
    btn_tm_set_v_Callback(hObject, eventdata, handles)
elseif flag_sweep == 1
    stopManualSpeedSelection(handles);
end
guidata(hObject, handles);


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


% --- Executes on button press in btn_mean_p.
function btn_mean_p_Callback(hObject, eventdata, handles)
flag_mean_p = getappdata(handles.output, 'flag_mean_p');
if flag_mean_p == 0
    flag_mean_p = 1;
    set(handles.btn_mean_p, 'BackgroundColor', 'green');
    datetime = clock;
    setappdata(handles.output, 'datetime_start_mean', datetime);
    setappdata(handles.output, 'datetime_start_mean_del', datetime);
    
    setappdata(handles.output, 'mean_p_n', 0);
    setappdata(handles.output, 'mean_p_p', 0);
elseif flag_mean_p == 1
    flag_mean_p = 0;
    set(handles.btn_mean_p, 'BackgroundColor', [0.94, 0.94, 0.94]);
    mean_p_p = getappdata(handles.output, 'mean_p_p');
    fprintf('mean p: %.3f [m]\n', mean_p_p);
end
setappdata(handles.output, 'flag_mean_p', flag_mean_p);
guidata(hObject, handles);
