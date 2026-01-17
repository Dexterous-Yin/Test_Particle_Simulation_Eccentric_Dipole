function B = Calc_B_ED(R, Phi, B_field, t, Scaling, params)
% calculates B-field at positions specified by array R
% B = [Bx; By; Bz] is a 2D array where first row is Bx-s, second row is
% By-s, third row is Bz-s.
% X = [x; y; z] is a 2D array where first row is x-s, second row is
% y-s, third row is z-s.
% This function operates on normalized
% quantities, i.e., x,y,z are in [R.E.], B in [m0*c^2/q/RE], t is in [R.E./c]
% =====================================================================
% This version is used to provide the magnetic field when tracing back
% the motion of particles. Because the dt-->-dt, so there is no need to
% reverse the direction of magnetic field.
% The ULF wave field is added the time factor
% created from Calc_B_3.m by yzf 20200318
% modified to use IGRF model. Zefan 20241219
% use igrfmex function. Zefan 20241225
% simplified for eccentric dipole. Zefan 20250327

Field_Model = B_field.Field_Model;
% BE = params.BE*1e9; % [nT] equatorial Earth B-field revised.

% x = R.*cos(Phi);
% y = R.*sin(Phi);
% z = zeros(size(R));
R_size = length(R);

if Field_Model == 0 % calculate B as analytical dipole
    % r = sqrt(x.^2 + y.^2 + z.^2);
    % r_3 = r.^3;
    % r_5 = r.^5;
    % 
    % % calculate field components in [B.E.]
    % B_x = - 3.*x.*z./r_5;
    % B_y = - 3.*y.*z./r_5;
    % B_z = 1./r_3 - 3*z.*z./r_5;
    % 
    % [B_phi,B_r] = cart2pol(B_x,B_y);
    % 
    % B = [B_r; B_phi; B_z];

    B_z = 1./R.^3;
    % [B_r; B_phi; B_z];
    B = [zeros(1,R_size);zeros(1,R_size);B_z];

else
    error('Calc_B:SettingError','Model Error.');
end

% Normalize B.
B = B/Scaling.B;
end

