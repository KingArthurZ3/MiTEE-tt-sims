function [tau_ctl] = calcControlTorque(controller_params,controller_inputs)
%calcControlTorque calculates the control torque produced given a set of
%parameters and inputs
%Inputs:
%   TBD
%Outputs:
%   tau_ctl - 1x3 vector representing control torque in x,y,z body frame

%              ------Function to be implemented here------

tau_ctl = [-0.01,24,60]*10^-3;
tau_ctl = [0 0 0];

end

