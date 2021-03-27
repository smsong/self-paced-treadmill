clear; clc;

i_gui = 2;
% 1: GUI_SPT
% 2: GUI_WalkTests

addpath('./fcn/');

dt_speedgoat=0.001; % 1000 hz

delete(timerfindall)

switch i_gui
    case 1
        h_gui = GUI_SPT();
    case 2
        h_gui = GUI_WalkTests();
end