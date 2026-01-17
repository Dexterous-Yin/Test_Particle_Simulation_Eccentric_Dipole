function params= set_parameters()
% fundamental parameters
params.RE = 6371.2e3; %[m]
% params.BE = 2.9744e4*1e-5; %[T]
params.ele_m = 9.1094e-31; %[kg]
params.ele_q = -1.6022e-19; %[C]
params.c = 2.99792e8; %[m/s]
params.keV = 1e3*abs(params.ele_q);

% drift and bounce functions
params.T0 = 1.3802;
params.T1 = 0.7405;




end



