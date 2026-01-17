%{ 
This script is used to trace particle's mapping point at equator in eccentric dipole field.
%}

% Basic settings
% eccentric dipole field settings
timenow = datetime(2024,04,09,21,00,00);
timenow_datenum = datenum(timenow);
[B0,thetan,phin,x0,y0,z0,trans_matrix] = recalc_ED_params(timenow_datenum);
% MEX geopack
igrfmex_wrapper_GC;

% PARAMETERS - DO NOT CHANGE - COUPLED TO CALCULATION PROCESS
params = set_parameters();
RE_km = params.RE/1e3;
% add eccentric params
params.BE = B0*1e-9; %[T]
params.x0 = x0;
params.y0 = y0;
params.z0 = z0;
params.trans_matrix = trans_matrix;

Earth_center_pos_mag = trans_matrix*[-x0,-y0,-z0]';
phi0 = atan2(Earth_center_pos_mag(1),Earth_center_pos_mag(2))/pi*180;
mirrphi = 90-phi0;

% spacecraft data
trange = ['2024-04-09/21:00:00';'2024-04-09/21:40:00'];
trange_datenum = datenum(trange,'yyyy-mm-dd/HH:MM:SS');
trange_datetime = datetime(trange_datenum,"ConvertFrom",'datenum');
pa_width = 10;
%
% load MSS data
load('PreCalc_A_0409.mat');

% particle tracing settings
E_set = energy_mid_MSS; %[keV]
pos_geo_set = pos_geo_in_Re;
sim_trange_set = time_date_MSS;
pa_set = 0:5:90; %[5:10:180];

% trans particle's position from GEO to approximate SM
L_set = zeros(size(sim_trange_set));
MLAT_set = zeros(size(sim_trange_set));
Phi_set = zeros(size(sim_trange_set));
Phi_mag_set = zeros(size(sim_trange_set));
PAeq_set = zeros(length(pa_set),length(sim_trange_set));
DLC_set = zeros(size(sim_trange_set));
DLC_Height = zeros(size(sim_trange_set));
BLC_set = zeros(size(sim_trange_set));
for Ti = 1:length(sim_trange_set)
    timestart = datetime(sim_trange_set(Ti),'ConvertFrom','datenum');
    iyear = year(timestart);
    imonth = month(timestart);
    iday = day(timestart);
    idoy = day(timestart,"dayofyear");
    ihour = hour(timestart);
    iminute = minute(timestart);
    isecond = second(timestart);
    recalc(iyear,idoy,ihour,iminute,isecond,-400,0,0);
    [sun_x_geo,sun_y_geo,sun_z_geo] = geogsw(1,0,0,-1);
    sun_pos_mag = trans_matrix*[sun_x_geo,sun_y_geo,sun_z_geo]';
    sun_pos_phi = atan2(sun_pos_mag(2),sun_pos_mag(1));
    % transformation from MAG to SM
    trans_matrix_sm = [cos(sun_pos_phi),sin(sun_pos_phi),0;...
        -sin(sun_pos_phi),cos(sun_pos_phi),0;...
        0,0,1];
    pos_ini_mag = trans_matrix*[pos_geo_set(1,Ti)-x0;pos_geo_set(2,Ti)-y0;pos_geo_set(3,Ti)-z0];
    pos_ini_sm = trans_matrix_sm*trans_matrix*[pos_geo_set(1,Ti)-x0;pos_geo_set(2,Ti)-y0;pos_geo_set(3,Ti)-z0];
    % get equatorial position based on dipole field
    R_ini_sm = sqrt(pos_ini_sm(1)^2+pos_ini_sm(2)^2+pos_ini_sm(3)^2);
    Lat_ini_sm = asin(pos_ini_sm(3)/R_ini_sm);
    L_ini_sm = R_ini_sm/cos(Lat_ini_sm)^2;
    Beq_ratio = cos(Lat_ini_sm)^6/sqrt(1+3*sin(Lat_ini_sm)^2);

    % record
    L_set(Ti) = L_ini_sm;
    Phi_set(Ti) = atan2(pos_ini_sm(2),pos_ini_sm(1));
    MLAT_set(Ti) = Lat_ini_sm/pi*180;
    
    Phi_mag_now = atan2(pos_ini_mag(2),pos_ini_mag(1));
    loss_slat = fzero(@(s) sqrt((L_ini_sm.*(1-s.^2).^(3/2).*cos(Phi_mag_now)-Earth_center_pos_mag(1)).^2+(L_ini_sm.*(1-s.^2).^(3/2).*sin(Phi_mag_now)-Earth_center_pos_mag(2)).^2+(L_ini_sm.*(1-s.^2).*s-Earth_center_pos_mag(3)).^2)-1-100/6371.2, 0);
    if loss_slat>0
        error('fzero error');
    end
    Bloss_ratio = sqrt(1+3*sin(Lat_ini_sm)^2)/cos(Lat_ini_sm)^6/(sqrt(1+3*loss_slat^2)/(1-loss_slat^2)^3);
    if Bloss_ratio>=1
        BLC_set(Ti) = 90;
    else
        BLC_set(Ti) = asin(sqrt(Bloss_ratio))/pi*180;
    end
    Phi_mag_set(Ti) = Phi_mag_now;

    for PAi=1:length(pa_set)
        PAeq_ini_sm = asin(sqrt(Beq_ratio*sind(pa_set(PAi))^2))/pi*180;
        PAeq_set(PAi,Ti) = PAeq_ini_sm;
    end

    % DLC
    L_now = L_ini_sm;
    if sqrt((L_now.*cosd(mirrphi)-Earth_center_pos_mag(1)).^2+(L_now.*sind(mirrphi)-Earth_center_pos_mag(2)).^2+(0-Earth_center_pos_mag(3)).^2)<=1+100/6371.2
        DLC_set(Ti) = 90;
    else
        sLat = fzero(@(s) sqrt((L_now.*(1-s.^2).^(3/2).*cosd(mirrphi)-Earth_center_pos_mag(1)).^2+(L_now.*(1-s.^2).^(3/2).*sind(mirrphi)-Earth_center_pos_mag(2)).^2+(L_now.*(1-s.^2).*s-Earth_center_pos_mag(3)).^2)-1-100/6371.2, 0);
        if ~isreal(sLat) || sLat>0
            error('check fzero');
        end
        sLat = abs(sLat);
        if abs(Lat_ini_sm)>=asin(sLat)
            DLC_set(Ti) = 90;
        else
            temp_ratio = (sqrt(1+3*sin(Lat_ini_sm)^2)/cos(Lat_ini_sm)^6)/(sqrt(1+3*sLat^2)/(1-sLat^2)^3);
            DLC_set(Ti) = asin(sqrt(temp_ratio))/pi*180;
        end
        DLC_Height(Ti) = (sqrt((L_now.*(1-sLat.^2).^(3/2).*cos(Phi_mag_now)-Earth_center_pos_mag(1)).^2+(L_now.*(1-sLat.^2).^(3/2).*sin(Phi_mag_now)-Earth_center_pos_mag(2)).^2+(L_now.*(1-sLat.^2).*(-sLat)-Earth_center_pos_mag(3)).^2)-1)*RE_km;
    end
end
MLT_set = Phi_set./pi*12+12;

%% Tracing part
% set tracing parameter sets
loopparams = [];
E_loop_range = 11; %1:length(energy_mid);
T_loop_range = 2011; %2011; %1771; %1891; %1:length(time_date_MSS);
PA_loop_range = 18; %1:length(pa_set);
for Ei = E_loop_range
    for PAi = PA_loop_range
        PAnow = pa_set(PAi);
        if PAnow>90
            PAnow = 180-PAnow;
        end
        inDLCid = find(BLC_set(T_loop_range)<PAnow);
        loopparams_now = [Ei+zeros(length(inDLCid),1),PAi+zeros(length(inDLCid),1),T_loop_range(inDLCid)'];
        loopparams = [loopparams;loopparams_now];
    end
end
loopparams_n = length(loopparams(:,1));

% initialization
Init_E = zeros(1,loopparams_n);
Init_L = zeros(1,loopparams_n);
Init_PA = zeros(1,loopparams_n);
Init_MLT = zeros(1,loopparams_n);
Min_Height = zeros(1,loopparams_n);

preview = 1;
height_bound = 100;
for loopi = 1:loopparams_n
    Ei = loopparams(loopi,1);
    PAi = loopparams(loopi,2);
    Ti = loopparams(loopi,3);

    % theoretical drift frequency
    theo_omega_d = omega_d(E_set(Ei)*params.keV,params.ele_m,params.ele_q,PAeq_set(PAi,Ti)/180*pi,L_set(Ti),params)+2*pi/24/60/60;
    % trace particle's properties, r, phi, PA
    [Particle,Fields,Time,Numerics] = Traj_setting(sim_trange_set,E_set,L_set,Phi_set,PAeq_set,Ei,PAi,Ti,theo_omega_d);
    % solve trajectory
    [Time_new, Res, Diagnostics] = Traj_particle_bounce_drift_ED(Particle, Fields, Time, Numerics, params);

    % output
    Time_old = Time;
    Time = Time_new;
    B_background = Diagnostics.B*params.BE*1e9; % [nT]
    B_norm = Diagnostics.norm_B*params.BE*1e9; %[nT]
    E_background = Diagnostics.E; %[mV/m]
    PA = Diagnostics.PA/pi*180;
    W_keV = Diagnostics.W*1e3;
    % SM coordinate
    R = Res(1,:);
    Phi = Res(2,:);
    MLT = Phi./pi*12+12;
    MLT = modify_MLT(MLT);
    Height = Diagnostics.Height*RE_km;
    GEOPOS = Diagnostics.GEOPOS;
    GEOPOSMIRR = Diagnostics.GEOPOSMIRR;
    GEOlon = atan2(GEOPOS(2,:),GEOPOS(1,:))./pi*180;
    GEOlonMIRR = atan2(squeeze(GEOPOSMIRR(2,2,:)),squeeze(GEOPOSMIRR(2,1,:)))./pi*180;
    GEOlat = asin(GEOPOS(3,:)/sqrt(GEOPOS(1,:).^2+GEOPOS(2,:).^2+GEOPOS(3,:).^2))./pi*180;
    Min_Height(loopi) = min(Height,[],'all');

    % store
    if isempty(Diagnostics.isloss) && min(Height,[],'all')>height_bound
        Init_E(loopi) = W_keV(end);
        Init_L(loopi) = R(end);
        Init_PA(loopi) = PA(end);
        Init_MLT(loopi) = MLT(end);
    end
    % if isempty(Diagnostics.isloss)
    %     break
    % end

    disp([num2str(loopi),'/',num2str(loopparams_n)]);
end

%% rearrange the results
InitialState.E = zeros(length(E_set),length(pa_set),length(sim_trange_set));
InitialState.L = zeros(length(E_set),length(pa_set),length(sim_trange_set));
InitialState.PA = zeros(length(E_set),length(pa_set),length(sim_trange_set));
InitialState.MLT = zeros(length(E_set),length(pa_set),length(sim_trange_set));

for loopi = 1:loopparams_n
    Ei = loopparams(loopi,1);
    PAi = loopparams(loopi,2);
    Ti = loopparams(loopi,3);
    InitialState.E(Ei,PAi,Ti) = Init_E(loopi);
    InitialState.L(Ei,PAi,Ti) = Init_L(loopi);
    InitialState.PA(Ei,PAi,Ti) = Init_PA(loopi);
    InitialState.MLT(Ei,PAi,Ti) = Init_MLT(loopi);
end

duration = 480;
E_PeakTime = datetime(2024,04,09,21,30,00)+seconds(duration/2)-seconds(540);
E_Amplitude = [-91.0000  320.0000 -275.2];
E_Parameters = [0.1,duration/2];  %[alpha,sigma]
E_ZeroTime = datetime(2024,04,09,20,00,00);
% save(['Result_PAD_',datestr(now,'mmddHHMM'),'.mat']);


%% PLOT PART
timestart = datetime(sim_trange_set(Ti),'ConvertFrom','datenum');
fig = figure('Color',[1 1 1]);
layout = tiledlayout(6, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
timestart_num = datenum(timestart);
Time_num = timestart_num+Time./60./60./24;
xlim_range = [datenum('2024-04-09/19:30:00'),timestart_num]; %[min(Time_num),timestart_num]; %[datenum('2024-04-09/20:00:00'),timestart_num];
if Ti == 1771
    xlim_range = [datenum('2024-04-09/20:25:00'),timestart_num]; 
end
if Ti == 2011
    xlim_range = [datenum('2024-04-09/20:05:00'),timestart_num]; 
end
xlim_id = find(Time_num>=xlim_range(1) & Time_num<=xlim_range(2));
xticks_datetime = dateshift(timestart+seconds(min(Time)),"start","hour","previous"):minutes(10):timestart;
xticks = datenum(xticks_datetime);
xticks_label = datestr(xticks,'HHMM');

% MAP
nexttile;
load coastlines
plot(coastlon,coastlat,'color',[0,0,0],'LineWidth',1); hold on;
% axis equal; 
grid on;
xlim([-180,180]); ylim([-60,60]);
% plot footpoint
Height_mirr = Height(:,xlim_id);
cbarmin = min(Height_mirr,[],"all");
cbarmax = max(Height_mirr,[],"all");
Einter = cbarmin:(cbarmax-cbarmin)/255:cbarmax;
cmap = colormap(turbo);
% mirror lat/lon
for southid = [2,1]
    south_mirr_pos = squeeze(GEOPOSMIRR(southid,:,:));
    south_mirr_lon = atan2(south_mirr_pos(2,:),south_mirr_pos(1,:))./pi*180;
    south_mirr_lat = asin(south_mirr_pos(3,:)./sqrt(south_mirr_pos(1,:).^2+south_mirr_pos(2,:).^2+south_mirr_pos(3,:).^2))./pi*180;
    south_mirr_height = Height_mirr(southid,:);
    south_mirr_id = xlim_id(1):300:xlim_id(end);
    if south_mirr_id(end) ~= xlim_id(end)
        south_mirr_id = [south_mirr_id,xlim_id(end)];
    end
    for linei = 1:length(south_mirr_id)-1
        tempcolor = interp1(Einter,cmap,mean(south_mirr_height(south_mirr_id(linei):south_mirr_id(linei+1))));
        % plot(south_mirr_lon(south_mirr_id(linei):south_mirr_id(linei+1)),south_mirr_lat(south_mirr_id(linei):south_mirr_id(linei+1)),'-','Color',tempcolor,'LineStyle','none','MarkerSize',6,'Marker','o','MarkerFaceColor',tempcolor,'MarkerEdgeColor',tempcolor); hold on;
        plot(south_mirr_lon(south_mirr_id(linei):south_mirr_id(linei+1)),south_mirr_lat(south_mirr_id(linei):south_mirr_id(linei+1)),'-','Color',tempcolor,'LineStyle','-','LineWidth',6); hold on;
    end
end
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',1.5,'FontSize', 14,'GridLineStyle','--');
% xlabel('Longitude (deg)'); 
ylabel('Latitude (deg)','FontSize',16);
cbar = colorbar('FontSize',14);
clim([cbarmin,cbarmax]);
cbar.Label.String = 'Height (km)'; %'Energy (keV)';
cbar.Label.FontSize = 16;
cbar.LineWidth = 1;
cbar.TickLength = 0.02;
set(gca,'XTick',[-180:30:150],'XTickLabelRotation',0,'YTick',[-60:20:60]);
title(string(timestart,'yyyy-MM-dd/HH:mm:ss'));
% set(gca,'TickDir','out','TickLength',[0.012,0.03]);
% set(gca,'XTickLabel',{});

% window control
% screen size
set(0, 'Units', 'pixels');
screenSize = get(0, 'ScreenSize');
screenHeight = screenSize(4);

% set the height of figure to the screen height
fig.Units = 'pixels'; 
fig.Position(2) = 0; 
fig.Position(4) = screenHeight;
fig.Position(3) = 480;
fig.Position(1) = (screenSize(3) - fig.Position(3)) / 2;

%
% PHI
nexttile
colororder({'k','b'})
yyaxis left
plot(Time_num,GEOlon,'k-','LineWidth',2);
set(gca,'YMinorTick','on'); ylabel('GEO Lon (deg)','FontSize',16);
if Ti == 2011
    ylim([-40,45]);
end
if Ti == 1771
    ylim([-40,25]);
end
yyaxis right
plot(Time_num,smoothdata(MLT,"movmean",1),'b-','LineWidth',2);
set(gca,'XTick',xticks,'XTickLabel',{});
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',1.5,'FontSize', 14, 'XTickLabel',{});
ylabel('MLT','FontSize',16); % ylim([0,24]); %ylim(refine(ylim));
ylim([14,24]);
if Ti == 1771 || Ti == 2011
    ylim([18,24]);
end
xlim(xlim_range);
% set(gca,'TickDir','out','TickLength',[0.012,0.03]);

% Height
nexttile
plot(xlim_range,[100,100],'k--','LineWidth',1.5); hold on;
% if Ti>0 %==1831
%     plot([min(Time_num),max(Time_num)],[height_A(T_loop_range),height_A(T_loop_range)],'k--','LineWidth',1);
%     text(Time_num(1), height_A(T_loop_range), [sprintf('%5.1f', height_A(T_loop_range)),'km'], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 14);
% end
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',1.5,'FontSize', 14, 'XTickLabel',{});
if Height(1,end)<=101
    plot(Time_num,smoothdata(Height(1,:),"movmean",1),'k-','LineWidth',2); hold on;
    ylabel('North Hmirr (km)','FontSize',16);
else
    plot(Time_num,smoothdata(Height(2,:),"movmean",1),'k-','LineWidth',2); hold on;
    ylabel('South Hmirr  (km)','FontSize',16);
    [minHeight,minid] = min(Height(2,:));
    if isempty(Diagnostics.isloss) 
        % plot([Time_num(minid),Time_num(minid)],[0,1000],'k--','LineWidth',1); hold on;
        text(Time_num(1), minHeight, [sprintf('%5.1f', minHeight),'km'], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 14);
    end
    
end
xlim(xlim_range);
if Ti == 1831
    ylim([50,600]);
    set(gca,'YTick',100:100:600,'YTickLabel',[100:100:600]);
end
if Ti == 1771
    ylim([50,500]);
end
if Ti == 2011
    ylim([50,700]);
end
set(gca,'XTick',xticks,'XTickLabel',{});
% set(gca,'TickDir','out','TickLength',[0.012,0.03]);

% Ephi
nexttile
plot(Time_num,smoothdata(E_background(2,:),"movmean",1),'k-','LineWidth',2); hold on;
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',1.5,'FontSize', 14);
ylabel('Ephi (mV/m)','FontSize',16);
xlim(xlim_range);
ylim([-5,0]);
if Ti == 2011
    ylim([-6.5,0]);
end
if Ti == 1771
    ylim([-3,0]);
end
set(gca,'XTick',xticks,'XTickLabel',{});
% set(gca,'TickDir','out','TickLength',[0.012,0.03]);

% W
nexttile
plot(Time_num,smoothdata(W_keV,"movmean",1),'k-','LineWidth',2); hold on;
% plot(Time(south_mirri),W_keV(south_mirri),'k-','LineWidth',2);
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',1.5,'FontSize', 14, 'XTickLabel',{});
ylabel('Energy (keV)','FontSize',16); % ylim(refine(ylim)); % ylim([99,101]); % ylim(refine(ylim));
xlim(xlim_range);
ylim([65,73]);
if Ti == 1771
    ylim([68,73]);
end
if Ti == 2011
    ylim([63,73]);
end
set(gca,'XTick',xticks,'XTickLabel',{});
% set(gca,'TickDir','out','TickLength',[0.012,0.03]);


% R
nexttile
plot(Time_num,smoothdata(R,"movmean",1),'k-','LineWidth',2);
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',1.5,'FontSize', 14);
ylabel('L','FontSize',16); 
ylim([1.86,1.96]);
if Ti == 1771
    ylim([1.95,1.99]);
end
if Ti == 2011
    ylim([1.74,1.88]);
end
xlim(xlim_range);
set(gca,'XTick',xticks,'XTickLabel',{});
set(gca,'XTick',xticks,'XTickLabel',xticks_label,'XTickLabelRotation',0);
% set(gca,'TickDir','out','TickLength',[0.012,0.03]);



%% FUNCTIONS
% drift omega
function res = omega_d(W,m,q,aeq,L,params)
gamma = 1+W./(m*params.c^2);
v = params.c*sqrt(1-1./gamma.^2);
y = sin(aeq);
T = params.T0-1/2*(params.T0-params.T1).*(y+y.^(1/2));
D = 1/12*(4*params.T0-(3*params.T0-5*params.T1).*y-(params.T0-params.T1).*(y.*log(y)+y.^(1/2)));
res = -3.*gamma.*m.*v.^2*L./(q.*params.BE.*params.RE.^2).*D./T;
end

function MLT = modify_MLT(MLT)
temp = find(MLT>24 | MLT<0);
while ~isempty(temp)
    MLT(MLT>24) = MLT(MLT>24)-24;
    MLT(MLT<0) = MLT(MLT<0)+24;
    temp = find(MLT>24 | MLT<0);
end
end

function [Particle,Fields,Time,Numerics] = Traj_setting(sim_trange_set,E_set,L_set,Phi_set,PAeq_set,Ei,PAi,Ti,theo_omega_d)

% initial condition
Particle.Charge = -1.0; %[elementary charge]
Particle.Mass = 1.0; %[el-n mass]
Particle.W_ini = E_set(Ei)/1e3;
Particle.R_ini = L_set(Ti);
Particle.Phi_ini = Phi_set(Ti);
Particle.Pitch_ini = PAeq_set(PAi,Ti);
Particle.y_ini = sind(Particle.Pitch_ini);

% initial time
timestart = datetime(sim_trange_set(Ti),'ConvertFrom','datenum');
iyear = year(timestart);
imonth = month(timestart);
iday = day(timestart);
idoy = day(timestart,"dayofyear");
ihour = hour(timestart);
iminute = minute(timestart);
isecond = second(timestart);

% parameters for electric field
duration = 480;
E_PeakTime = datetime(2024,04,09,21,29,00)-seconds(duration/2); %datetime(2024,04,09,21,30,00)+seconds(duration/2)-seconds(540);
E_Amplitude = [-91.0000  320.0000 -275.2]; %[mV/m]
E_Parameters = [0.1,duration/2];  %[alpha,sigma]
E_ZeroTime = datetime(2024,04,09,20,00,00); %E_PeakTime-seconds(duration/2)-seconds(100); %E_PeakTime-hours(1); %datetime(2024,04,09,21,30,00)-seconds(duration)-seconds(100); %E_PeakTime-hours(1);

% background fields
Fields.B_field.Field_Model  = 0; % 0 - analytical dipole, 1 - IGRF
Fields.E_field.Time_Peak  = E_PeakTime; % variable for E field setup
Fields.E_field.Time_SecondLeap = (Fields.E_field.Time_Peak-timestart)/seconds(1);
Fields.E_field.Amp_Peak  = E_Amplitude; % [mV/m]
Fields.E_field.Parameters = E_Parameters;
Fields.E_field.Time_Start  = timestart;
Fields.E_field.IYEAR_Start  = iyear;
Fields.E_field.IDOY_Start  = idoy;
Fields.E_field.IHOUR_Start  = ihour;
Fields.E_field.IMINUTE_Start  = iminute;
Fields.E_field.ISECOND_Start  = isecond;
Fields.E_field.V_con = 1; % nonzero: open; 0: off

Numerics.X_AbsTol = 1.e-12; % [RE]

duration = (datetime(sim_trange_set(Ti),'ConvertFrom','datenum')-E_ZeroTime)/seconds(1)+3*pi/theo_omega_d;

Time = 0:-0.5:-duration; % output particle position at specified moments of time
end