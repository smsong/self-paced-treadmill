function [v_tm_tgt, a_tm_tgt] = speed_controller(...
    flag_step, com_step, dcom_step, v_tm, ...
    p_off, G_p, G_v, r_v_tm_tgt, del_t_tgt, ...
    v_tm_tgt0, a_tm_tgt0)
% PI controller

a_min = 0.001; % minimum acceleration

if ~flag_step
    v_tm_tgt = v_tm_tgt0;
    a_tm_tgt = a_tm_tgt0;
    return;
end

del_v_tgt = G_p*(com_step-p_off) + G_v*dcom_step;
a_tgt = abs(del_v_tgt/del_t_tgt);

v_tm_tgt = v_tm + del_v_tgt;

v_tm_tgt = min(v_tm_tgt, r_v_tm_tgt(2));
v_tm_tgt = max(v_tm_tgt, r_v_tm_tgt(1));
a_tm_tgt = max(a_tgt, a_min);