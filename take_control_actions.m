function ps = take_control_actions(ps,sub_grids,ramp_rate,dt,it_no,opt)
% usage: ps = take_control_actions(ps,sub_grids,ramp_rate,dt,it_no,opt)
%
% Compute and implement emergency control actions.
%  Interfaces with comm model, if requested in the options

switch opt.sim.control_method
    case 'emergency_control'
        ps = old_control_actions(ps,sub_grids,ramp_rate,dt,it_no,opt);
    case 'mpc'
        % call Pooya's code
    case 'distributed_mpc'
        
    case 'none'
        % Do nothing
        %ps = ps;
    otherwise 
        error('Undefined control method');
end
