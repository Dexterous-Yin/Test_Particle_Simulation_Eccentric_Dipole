function E = Calc_E_ED(R, Phi, E_field, t, Scaling, params)

% created by yzf 20200318
% This version is used to provide the electeic field when tracing back
% the motion of particles, with the growing factor for ULF wave field.
% MODIFICATION
% used in igrf model. first add corotation electric field. Zefan 20241225
% add convection electric field

Time_Start = E_field.Time_Start;
V_con = E_field.V_con;

PeakTime = E_field.Time_Peak;
PeakTime_Second = E_field.Time_SecondLeap;
PeakAmp = E_field.Amp_Peak;
EParameters = E_field.Parameters;

RE = params.RE; %[m]
BE = params.BE; %[T]

R_size = size(R,2);
T_size = size(t,2); % should be 1 or the same size of position

E_background=zeros(3,R_size);
if V_con ~=0
    simtime_leap = t*Scaling.Time; %[s]
    % E_phi
    E_background(2,:) = calc_E_phi(R,Phi,simtime_leap,PeakTime_Second,PeakAmp,EParameters);
end

E=E_background;

% Normalize E.
E = E/Scaling.E;
end

function E_phi = calc_E_phi(r,phi,timenow,timepeak,amppeak,Eparams)
tempt = timenow-timepeak;
alpha = Eparams(1);
sigma = Eparams(2);
a = amppeak(1);
b = amppeak(2);
c = amppeak(3);
rbound_1 = (-b-sqrt(b^2-4*a*c))/(2*a);
rbound_2 = (-b+sqrt(b^2-4*a*c))/(2*a);
rbound_min = min([rbound_1,rbound_2]);
rbound_max = max([rbound_1,rbound_2]);
E_phi = (a.*r.^2+b.*r+c).*cos(phi).*(tanh(alpha.*(tempt+sigma))+tanh(alpha.*(-tempt+sigma)))./(2.*tanh(alpha*sigma));
E_phi(r>rbound_max) = 0;
E_phi(r<rbound_min) = 0;
end