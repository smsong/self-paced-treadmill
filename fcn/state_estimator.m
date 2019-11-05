function [flag_step, com_step, dcom_step, memory] = ...
    state_estimator(flag_reset, time, grfs_raw, v_tm, memory)
% Kalman Filter

grfs_queue_n = 10; % length of sensor queue to filter grf data

i = 1;
time_0 = memory(i); i = i + 1; % i == 1
x = memory(i+(0:1)); i = i + 2; % i == 2, 3
P = reshape(memory(i+(0:3)), [2, 2]); i = i + 4; % i == 4, 5, 6, 7
stance_L_0 = memory(i); i = i + 1; % i == 8
stance_R_0 = memory(i); i = i + 1; % i == 9
len_step = memory(i); i = i + 1;
leg_step_1 = memory(i); i = i + 1; % i == 11
leg_step_0 = memory(i); i = i + 1;
p_step_1 = memory(i); i = i + 1; % i == 13
p_step_0 = memory(i); i = i + 1; % i == 14
t_step_1 = memory(i); i = i + 1;
t_step_0 = memory(i); i = i + 1;
n_step = memory(i); i = i + 1; % i == 17
com_step = memory(i); i = i + 1;
dcom_step = memory(i); i = i + 1;
v_tm_mean = memory(i); i = i + 1;
n_ts = memory(i); i = i + 1;
grfs_queue_i = memory(i); i = i + 1;
grfs_queue_buffer = reshape(memory(i+(0:grfs_queue_n*6-1)), [grfs_queue_n, 6]); i = i + grfs_queue_n*6;
% i - 1 == 32

%=======%
% RESET %
%=======%
if flag_reset
    time_0 = time;
    x = [0; 0]; % state: com, dcom
    P = [.0035 .0015; .0015 .0016]; % a posteriori error covariance matrix (from experiment (mocap as reference data))

    stance_L_0 = 0; % left stance
    stance_R_0 = 0; % right stance
    len_step = .7; % step length
    leg_step_1 = 0; % the most recent step was (0: L; 1: R)
    leg_step_0 = 1; % the step before that was (0: L; 1: R)
    p_step_1 = 0; % y position of the most recent step
    p_step_0 = 0; % y position of the previous step
    t_step_1 = 0; % time when the most recent step was made
    t_step_0 = -.5; % time when the previous step was made
    n_step = 0; % total number of steps
    com_step = 0; % mean com during step
    dcom_step = 0; % mean dcom during step
    v_tm_mean = 0; % mean v_tm during step
    n_ts = 0; % number of time steps since last step

    % sensor queue to filter grf data
    grfs_queue_i = 0;
    grfs_queue_buffer = zeros(grfs_queue_n, 6);
end

del_t = time - time_0;


%=================%
% FILTER GRF DATA %
%=================%
grfs_queue_i = grfs_queue_i + 1;
if grfs_queue_i > grfs_queue_n
    grfs_queue_i = grfs_queue_i - grfs_queue_n;
end
grfs_queue_buffer(grfs_queue_i,:) = grfs_raw;
grfs_filt = mean(grfs_queue_buffer);


%===============================%
% [KALMAN FILTER] PREDICT STATE %
%===============================%
u = [grfs_filt(1) + grfs_filt(4)]; % input (acceleration)
A = [1 del_t; 0 1];
B = [(del_t^2)/2; del_t];
noise_q2 = 2.9; % square of q: from experiment (mocap as reference data)
Q = B*B'*noise_q2;

x = A*x + B*u;
P = A*P*A' + Q;


%==================================%
% [KALMAN FILTER] NEW MEASUREMENT? %
%                 i.e. new step?   %
%==================================%
h_fp = 0.015; % forceplate height [m]
off_trdml = -.5*1.8034; % treadmill offset [m] (plate length = 71 inch = 1.8034 m)
th_grfz = .3*9.81; % grfz threshold for step detection [m/s^2] grfz = 9.81 -> full body

fyl_filt = grfs_filt(1);
fzl_filt = grfs_filt(2);
mxl_filt = grfs_filt(3);
fyr_filt = grfs_filt(4);
fzr_filt = grfs_filt(5);
mxr_filt = grfs_filt(6);

% update mean values of com, dcom, v_tm
com_step = (com_step*n_ts + x(1))/(n_ts + 1);
dcom_step = (dcom_step*n_ts + x(2))/(n_ts + 1);
v_tm_mean = (v_tm_mean*n_ts + v_tm)/(n_ts + 1);
n_ts = n_ts + 1;

% detect new step
flag_step = 0;
if stance_L_0 == 0 && fzl_filt > th_grfz % left step?
    flag_step = 1;
    p_step_0 = p_step_1;
    p_step_1 = -(h_fp*fyl_filt + mxl_filt)/fzl_filt + off_trdml;
    leg_step_1 = 0; % 0: L
end
if stance_R_0 == 0 && fzr_filt > th_grfz % right step?
    flag_step = 1;
    p_step_0 = p_step_1;
    p_step_1 = -(h_fp*fyr_filt + mxr_filt)/fzr_filt + off_trdml;
    leg_step_1 = 1; % 1: R
end

% if new step made
if flag_step
    t_step_0 = t_step_1;
    t_step_1 = time;
    len_step = p_step_1 - p_step_0;
    leg_step_0 = leg_step_1;
    n_step = n_step + 1;
    n_ts = 0;
end

% update stance data
stance_L_0 = fzl_filt > th_grfz;
stance_R_0 = fzr_filt > th_grfz;
p_step_0 = p_step_0 - v_tm*del_t;
p_step_1 = p_step_1 - v_tm*del_t;


%==================================%
% [KALMAN FILTER] UPDATE STATE
%                 if new step made %
%==================================%
% do not update when step is suspicious (i.e. stepping on same belt)
if flag_step && ((t_step_1 - t_step_0) < 1.2)
    flag_step_ok = 1;
else
    flag_step_ok = 0;
end
if flag_step_ok
    % measure output (y)
    com_mes = (p_step_1 + p_step_0)/2;
    dcom_mes = len_step/(t_step_1 - t_step_0) - v_tm_mean;
    y_mes = [com_mes; dcom_mes];

    % update state
%     C = [1 0; 0 1]; D = [0; 0];
    noise_v = [.025 0; 0 .085]; % from experiment (mocap as reference data)
    R = noise_v*noise_v';

    y_est = [com_step; dcom_step]; % y_est = C*[x(1); dcom_step] + D*u; % estimate output

    K = P/(P + R); % K = P*C'/(C*P*C' + R);
    x = x + K*(y_mes - y_est);
    P = (eye(2) - K)*P; % P = (eye(2) - K*C)*P;
else
    y_mes = [-100; -100];
end

u_y = [u; y_mes];

%==============================%
% SAVE DATA FOR NEXT ITERATION %
%==============================%
i = 1;
memory(i)= time; i = i + 1;
memory(i+(0:1)) = x; i = i + 2;
memory(i+(0:3)) = P(:); i = i + 4;
memory(i) = stance_L_0; i = i + 1; disp(memory(8))
memory(i) = stance_R_0; i = i + 1;
memory(i) = len_step; i = i + 1;
memory(i) = leg_step_1; i = i + 1;
memory(i) = leg_step_0; i = i + 1;
memory(i) = p_step_1; i = i + 1;
memory(i) = p_step_0; i = i + 1;
memory(i) = t_step_1; i = i + 1;
memory(i) = t_step_0; i = i + 1;
memory(i) = n_step; i = i + 1;
memory(i) = com_step; i = i + 1;
memory(i) = dcom_step; i = i + 1;
memory(i) = v_tm_mean; i = i + 1;
memory(i) = n_ts; i = i + 1;
memory(i) = grfs_queue_i; i = i + 1;
memory(i+(0:grfs_queue_n*6-1)) = grfs_queue_buffer; i = i + grfs_queue_n*6;
