function ps = take_control_actions(ps,sub_grids,ramp_rate,dt,it_no,opt)
% usage: ps = take_control_actions(ps,sub_grids,ramp_rate,dt,it_no,opt)
%
% Compute and implement emergency control actions.
%  Interfaces with comm model, if requested in the options

switch opt.sim.control_method
    case 'emergency_control'
        ps = old_control_actions(ps,sub_grids,ramp_rate,dt,it_no,opt);
    case 'emergency_control_dec' % decentralized emergency control (no MPC)
        ramp_limits = ramp_rate * dt;
        ps = dec_control(ps,sub_grids,ramp_limits,opt);
    case 'distributed_mpc'
        ramp_limits = ramp_rate * dt;
        ps = dist_mpc_control(ps,sub_grids,ramp_limits,opt);
    case 'none'
        % Do nothing
    otherwise 
        error('Undefined control method');
end
