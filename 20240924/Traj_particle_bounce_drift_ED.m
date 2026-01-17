function [Time_new, Res, Diagnostics] = Traj_particle_bounce_drift_ED(Particle, Fields, Time, Numerics, params)
%  Calculates particle's trajectory at times specified by Time[]
%  array.
%  Just focus on bounce drift motion.
%  Based on eccentric dipole magnetic field. 250326 Zefan.
% Returns:
%      1st Res(1,:) = R, [RE] - L-shell
%      2nd Res(2,:) = phi, azimuthal angle
%      3rd Res(3,:) = y, sin(PAeq)

% caluclate normalization coefficients
Scaling = Calc_Scaling(Particle, params);

% normalize input parameters
[Particle, Fields, Time] = Normalize(Particle, Fields, Time, Scaling);

% calculate initial magnetic momentum (invariant)
W0 = Particle.W_ini;
B0 = norm(Calc_B_ED(Particle.R_ini, Particle.Phi_ini, Fields.B_field, Time(1), Scaling, params));
Pitch0 = Particle.Pitch_ini;
M0 = W0./B0.*(W0/2 + 1).*(sin(Pitch0)).^2 ; % magnetic momentum

X_ini = [Particle.R_ini; Particle.Phi_ini; Particle.y_ini]; % inital condition 4-vector [x;y;z;p_par]

% set particle position tolerances
AbsTol(1) = Numerics.X_AbsTol; % abs.tol. in R
AbsTol(2) = Numerics.X_AbsTol; % abs.tol. in Phi
AbsTol(3) = Numerics.X_AbsTol; % abs.tol. in y
% options = odeset('RelTol',2.5e-14,'AbsTol',AbsTol,'OutputFcn',@odetpbar);
% options = odeset('RelTol',2.5e-14,'AbsTol',AbsTol,'OutputFcn',@odetpbar,'Events',@(t,y) events(t,y,Fields,Scaling,params));
options = odeset('MaxStep', 500, 'RelTol',2.23e-14,'AbsTol',AbsTol,'Events',@(t,y) events(t,y,Fields,Scaling,params));
% options = odeset('RelTol',2.23e-14,'AbsTol',AbsTol,'Events',@(t,y) events(t,y,Fields,Scaling,params));


% CALCULATATION PART
Charge = Particle.Charge;
[t, X, te, Xe, ie] = ode45(@(t,X) Calc_dX_dt(t, X, M0, Charge, Fields, Scaling, params), Time, X_ini, options);

% scale variables back to those with units
X = X'; %transpose solution array
% calculate outputs for diagnostics
[Diagnostics.B, Diagnostics.norm_B, Diagnostics.E, Diagnostics.W, Diagnostics.gamma, ...
    Diagnostics.mu, Diagnostics.p_perp, Diagnostics.p_par, Diagnostics.PA, Diagnostics.Height, Diagnostics.GEOPOS, Diagnostics.GEOPOSMIRR] = Calc_Diagnostics(X, M0, Fields, t', Scaling, params); % calculate outputs
% calculate kinetic energy array
Diagnostics.W = Diagnostics.W*Scaling.Energy;
% scale parallel momentum back to w/units
Diagnostics.p_perp = Diagnostics.p_perp*Scaling.p_par;
Diagnostics.p_par = Diagnostics.p_par*Scaling.p_par;
% magnetic field along the trajectory
Diagnostics.B = Diagnostics.B*Scaling.B;
Diagnostics.norm_B = Diagnostics.norm_B*Scaling.B;
% electric field along the trajectory
Diagnostics.E = Diagnostics.E*Scaling.E;
% record the terminal contion
Diagnostics.isloss = te;

% check the boundary here.
if ~isempty(te)
    disp(['The boundary has been hit at time: ', num2str(te*Scaling.Time), ' s.']);
end

% scale position vectors back to w/units
Res(1,:) = X(1,:)*Scaling.r;
Res(2,:) = X(2,:);

Time_new = t'*Scaling.Time;

function dX_dt = Calc_dX_dt(t, X, M, Charge, Fields, Scaling, params)
% this is the function used by the ODE solver
% calculates r.h.s. of the g.c. drift equations
% X = [R;Phi;y], dX_dt = [dR/dt; dPhi/dt; dy/dt]
% operates on normalized quantities
R = X(1,1);
Phi = X(2,1);
y = X(3,1);
q = sign(Charge);
% magnetic field
B = Calc_B_ED(R, Phi, Fields.B_field, t, Scaling, params);
norm_B = norm(B); % magnitude of B
p_perp_sq = 2*M*norm_B;
p_par_sq = p_perp_sq*(1-y^2)/y^2;
gamma = sqrt(p_perp_sq + p_par_sq + 1); % gamma
B_2 = norm_B^2; % magnitude of B squared
% electric field
E = Calc_E_ED(R, Phi, Fields.E_field, t, Scaling, params); % E-field vector
% drift related functions
[u_m,Y,T,D] = omega_d(gamma,y,R,params);
u_e = cross(E,B)/B_2;

dR_dt = u_e(1);
dPhi_dt = u_e(2)/R + u_m/Scaling.Frequency + 2*pi/24/60/60/Scaling.Frequency;
dy_dt = -Y/4/T*y/R*dR_dt;

% build return vector dX_dt
dX_dt = [dR_dt; dPhi_dt; dy_dt];

function Scaling = Calc_Scaling(Particle,params)
% calculates scaling coefficients for the variables
B0 = params.BE*1e4; % [G]  Earth's B-field at the surface (equator)
RE = 6.3712e8; % [cm] - Earth's radius
m0 = 9.1094e-28; % [g] - electron mass
e  = 4.8032e-10; % [statcoul] - elementary charge
c  = 2.99792e10; % [cm/sec] - speed of light
MeV = 1.6022e-6; % [erg] - energy [erg] equivalent to 1 MeV
mV_m = 3.33564e-8; % [statvolt/cm] - E-field [statV/cm] equavalent to 1 mV/m
Mass = Particle.Mass;
Charge = abs(Particle.Charge);

Scaling.B = Mass*m0*c^2/Charge/e/RE/B0; % [B.E.]
Scaling.r = 1.0;    % [R.E.]
Scaling.Energy = Mass*m0*c^2/MeV; % [MeV]
Scaling.E = Mass*m0*c^2/Charge/e/RE/mV_m;   % [mV/m]
Scaling.Frequency = c/RE; % [Hz]
Scaling.Time = RE/c;    % [sec]
Scaling.Angle = 180/pi; % [deg/rad]
Scaling.p_par = 1; %[m0*c]

function [Particle_out, Fields_out, Time_out] = Normalize(Particle, Fields, Time, Scaling)
% normalizes the quantities by dividing them by Scaling coefficients
% Scaling - structure, containing scaling coefficients, generetaed by the 'Scaling fucntion'

Particle.W_ini = Particle.W_ini/Scaling.Energy;
Particle.Pitch_ini = Particle.Pitch_ini/Scaling.Angle;
Particle.R_ini = Particle.R_ini/Scaling.r;

Time = Time/Scaling.Time;

Particle_out = Particle;
Fields_out = Fields;
Time_out = Time;

function [B, norm_B, E, W, gamma, Mu, p_perp, p_par, PA, Height, GEOPOS, GEOPOSMIRR] = Calc_Diagnostics(X, M, Fields, Time, Scaling, params)
% calculates kinetic energy
% operates on normalized quantities
R = X(1,:);
Phi = X(2,:);
y = X(3,:);
B = Calc_B_ED(R,Phi, Fields.B_field, Time, Scaling, params);
E = Calc_E_ED(R,Phi, Fields.E_field, Time, Scaling, params);
norm_B = sqrt(B(1,:).^2 + B(2,:).^2 + B(3,:).^2);
p_perp = sqrt(2.0*M*norm_B);
PA = asin(y);
p_par = p_perp./tan(PA);
W = sqrt(1 + p_perp.^2 + p_par.^2) - 1.0;
gamma = W+1;
Mu = p_perp.^2./2./norm_B;
% get mirror height and GEO position
Height = zeros(2,length(R));
GEOPOS = zeros(3,length(R));
GEOPOSMIRR = zeros(2,3,length(R));
for Ti = 1:length(R)
    L_now = R(Ti);
    Phi_now = Phi(Ti);
    [Earth_center_sm,trans_matrix_sm] = get_Earth_center_pos(Time(Ti),Fields.E_field,Scaling,params);
    % mirror height
    [sLat,fval,exitflag,output] = fzero(@(s) sqrt(1+3*s^2)/(1-s^2)^3-1/y(Ti)^2, [0,0.99]);
    sLat = abs(sLat);
    cLat = sqrt(1-sLat^2);
    pos_mirr = [L_now*cLat^3*cos(Phi_now),L_now*cLat^3*sin(Phi_now),L_now*cLat^2*sLat; ...
        L_now*cLat^3*cos(Phi_now),L_now*cLat^3*sin(Phi_now),-L_now*cLat^2*sLat];
    height_mirr = ((pos_mirr(:,1)-Earth_center_sm(1)).^2 + (pos_mirr(:,2)-Earth_center_sm(2)).^2 + (pos_mirr(:,3)-Earth_center_sm(3)).^2).^0.5-1;
    Height(:,Ti) = height_mirr;
    % GEO position
    pos_now = R(Ti).*[cos(Phi(Ti));sin(Phi(Ti));0];
    pos_mag = trans_matrix_sm'*pos_now;
    pos_geo_ecc = params.trans_matrix'*pos_mag;
    pos_geo = pos_geo_ecc+[params.x0;params.y0;params.z0];
    GEOPOS(:,Ti) = pos_geo;
    % GEO position of mirror point
    for i=1:2
        posmirr_mag = trans_matrix_sm'*pos_mirr(i,:)';
        posmirr_geo_ecc = params.trans_matrix'*posmirr_mag;
        posmirr_geo = posmirr_geo_ecc+[params.x0;params.y0;params.z0];
        GEOPOSMIRR(i,:,Ti) = posmirr_geo;
    end
    
end

function [omega_d,Y,T,D] = omega_d(gamma,y,L,params)
v = params.c*sqrt(1-1./gamma.^2);
Y = 2*(1-y)*params.T0+(params.T0-params.T1)*(y.*log(y)+2*y-2*y.^(1/2));
T = params.T0-1/2*(params.T0-params.T1).*(y+y.^(1/2));
D = 1/12*(4*params.T0-(3*params.T0-5*params.T1).*y-(params.T0-params.T1).*(y.*log(y)+y.^(1/2)));
omega_d = -3.*gamma.*params.ele_m.*v.^2*L./(params.ele_q.*params.BE.*params.RE.^2).*D./T;

function [value,isterminal,direction] = events(t, y, Fields, Scaling, params)
% Locate the boundary the particle crosses and stop integration.
R_bound = 1+100e3/params.RE; % 100km altitude
% get the position of Earth's center
[Earth_center_sm,~] = get_Earth_center_pos(t,Fields.E_field,Scaling,params);
% get mirror point
Phi = y(2,:);
L = y(1,:);
sLat = fzero(@(s) sqrt(1+3*s^2)/(1-s^2)^3-1/y(3)^2, [0,0.99]);
cLat = sqrt(1-sLat^2);
pos_mirr = [L*cLat^3*cos(Phi),L*cLat^3*sin(Phi),L*cLat^2*sLat; ...
    L*cLat^3*cos(Phi),L*cLat^3*sin(Phi),-L*cLat^2*sLat];
value_mirr = ((pos_mirr(:,1)-Earth_center_sm(1)).^2 + (pos_mirr(:,2)-Earth_center_sm(2)).^2 + (pos_mirr(:,3)-Earth_center_sm(3)).^2).^0.5 - R_bound;
value = min(value_mirr);
if ~isreal(value)
    error('fzero failed.');
end
isterminal = 1;   % stop the integration
direction = 0;   % Detect all zero crossings
% direction = -1;   % negative direction
% direction = 1;   % positive direction

function [pos_sm,trans_matrix_sm] = get_Earth_center_pos(t,E_field,Scaling,params)
iyear = E_field.IYEAR_Start;
idoy = E_field.IDOY_Start;
ihour = E_field.IHOUR_Start;
iminute = E_field.IMINUTE_Start;
isecond = E_field.ISECOND_Start;
simtime_leap = t*Scaling.Time; %[s]
isecond_now = isecond+simtime_leap;
% recalc geopack parameters
igrfmex(1,iyear,idoy,ihour,iminute,isecond_now,-400,0,0);
[sun_x_geo,sun_y_geo,sun_z_geo] = igrfmex(8,1,0,0,-1);
sun_pos_mag = params.trans_matrix*[sun_x_geo,sun_y_geo,sun_z_geo]';
sun_pos_phi = atan2(sun_pos_mag(2),sun_pos_mag(1));
% transformation from MAG to SM
trans_matrix_sm = [cos(sun_pos_phi),sin(sun_pos_phi),0;...
    -sin(sun_pos_phi),cos(sun_pos_phi),0;...
    0,0,1];
pos_sm = trans_matrix_sm*params.trans_matrix*[-params.x0;-params.y0;-params.z0];