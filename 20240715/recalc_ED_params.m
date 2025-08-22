% This function is used to calcualte the IGRF params at a given time, and
% then estimate the parameters for the eccentric dipole magnetic field.
function [B0,thetan,phin,x0,y0,z0,trans_matrix] = recalc_ED_params(Time_num)
% IGRF parameters
% g10,g11,h11,g20,g21,h21,g22,h22
% IGRF_2000 = [-29619.4,-1728.2,5186.1,-2267.7,3068.4,-2481.6,1670.9,-458.0];
IGRF_2020 = [-29403.41,-1451.37,4653.35,-2499.78,2981.96,-2991.72,1676.85,-734.62];
IGRF_2025 = [-29350.0,-1410.3,4545.5,-2556.2,2950.9,-3133.6,1648.7,-814.2];
F2 = (Time_num-datenum('2020-01-01/00:00:00'))/365.25/5;
F1 = 1.-F2;
IGRF_params = IGRF_2020*F1+IGRF_2025*F2;
IGRF_params_cell = num2cell(IGRF_params);
[g10,g11,h11,g20,g21,h21,g22,h22] = IGRF_params_cell{:};

% calculate eccentric dipole
% dipole tilt
B0 = sqrt(g10^2+g11^2+h11^2);
c11 = sqrt(g11^2+h11^2);
cos_thetan = -g10/B0;
sin_thetan = c11/B0;
thetan = atan2(sin_thetan,cos_thetan);
cos_phin = -g11/c11;
sin_phin = -h11/c11;
phin = atan2(sin_phin,cos_phin);

% eccentric coordinates
L0 = 2*g10*g20+sqrt(3)*(g11*g21+h11*h21);
L1 = -g11*g20+sqrt(3)*(g10*g21+g11*g22+h11*h22);
L2 = -h11*g20+sqrt(3)*(g10*h21-h11*g22+g11*h22);
E = (L0*g10+L1*g11+L2*h11)/(4*B0^2);
x0 = (L1-g11*E)/(3*B0^2);
y0 = (L2-h11*E)/(3*B0^2);
z0 = (L0-g10*E)/(3*B0^2);

% transformation from eccentric GEO to MAG
trans_matrix = [cos_thetan*cos_phin,cos_thetan*sin_phin,-sin_thetan;...
    -sin_phin,cos_phin,0;...
    sin_thetan*cos_phin,sin_thetan*sin_phin,cos_thetan];

% eccentricity
x0 = -0.092739;
y0 = 0.086886;
z0 = 0.096215;
end