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

Eparams = params.Eparams;

RE = params.RE; %[m]
BE = params.BE; %[T]

R_size = size(R,2);
T_size = size(t,2); % should be 1 or the same size of position

E_background=zeros(3,R_size);
if V_con ~=0
    simtime_leap = t*Scaling.Time; %[s]
    % E_phi
    % [E_phi,E_r] = calc_E_y(R,Phi,simtime_leap,PeakTime_Second,PeakAmp,EParameters);
    E_phi = calc_E_phi(R,Phi,simtime_leap,PeakTime_Second,PeakAmp,EParameters);
    E_r = zeros(size(E_phi));

    % ULF wave
    picktimenum = datenum(Time_Start)+simtime_leap/60/60/24;
    ULF_E_Phi = calc_ULFE(R,Phi,picktimenum,Eparams);
    E_phi = E_phi+ULF_E_Phi;

    E_background(1,:) = E_r;
    E_background(2,:) = E_phi;
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

% function [E_phi,E_r] = calc_E_y(r,phi,timenow,timepeak,amppeak,Eparams)
% tempt = timenow-timepeak;
% alpha = Eparams(1);
% sigma = Eparams(2);
% k = amppeak(2);
% E0 = amppeak(1)+k.*(r-1);
% Ey = E0.*(tanh(alpha.*(tempt+sigma))+tanh(alpha.*(-tempt+sigma)))./(2.*tanh(alpha*sigma));
% E_phi = Ey.*cos(phi);
% E_r = zeros(size(r)); %Ey.*sin(phi);
% end

% function Ephi = calc_ULFE(L,Phi,Time,Eparams)
% Ephi = Eparams.E0.*sin(Eparams.ULF_omega.*(Time-Eparams.t0)*24*60*60-Eparams.ULF_m.*Phi+Eparams.phi0);
% end


% function E_phi = calc_E_phi(r,phi,timenow,timepeak,amppeak,Eparams)
% tempt = timenow-timepeak;
% sigma = Eparams(1);
% E_phi = amppeak.*exp(-tempt.^2/2/sigma^2).*cos(phi);
% end

function Ephi = calc_ULFE(L,Phi,Time,Eparams)
Ephi = zeros(size(Time));
preid = find(Time-Eparams.peaktime<0);
Ephi(preid) = Eparams.E0.*exp(-((Time(preid)-Eparams.peaktime)*24*60*60).^2./Eparams.tau1.^2).*cos(Eparams.ULF_m.*(Phi(preid)-Eparams.Phi0)-Eparams.ULF_omega.*(Time(preid)-Eparams.t0)*24*60*60+Eparams.phase0);
postid = find(Time-Eparams.peaktime>=0);
Ephi(postid) = Eparams.E0.*exp(-((Time(postid)-Eparams.peaktime)*24*60*60).^2./Eparams.tau2.^2).*cos(Eparams.ULF_m.*(Phi(postid)-Eparams.Phi0)-Eparams.ULF_omega.*(Time(postid)-Eparams.t0)*24*60*60+Eparams.phase0);
end