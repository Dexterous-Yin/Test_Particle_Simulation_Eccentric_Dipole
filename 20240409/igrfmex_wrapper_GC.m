%{
This is a MATLAB wapper for the igrfmex functions.
#define T04_S_IDX 0
#define RECALC_08_IDX 1
#define IGRF_GSW_08_IDX 2
#define SMGSW_08_IDX 3
#define MAGSM_08_IDX 4
#define GEIGEO_08_IDX 5
#define GEOMAG_08_IDX 6
#define GSWGSE_08_IDX 7
#define GEOGSW_08_IDX 8
%}

t04 = @(parmod,ps,x,y,z)igrfmex(0,parmod,ps,x,y,z);

% SUBROUTINE RECALC_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,VGSEX,VGSEY,VGSEZ)
% PREPARES ELEMENTS OF ROTATION MATRICES FOR TRANSFORMATIONS OF VECTORS BETWEEN SEVERAL COORDINATE SYSTEMS
% PREPARES COEFFICIENTS USED IN THE CALCULATION OF THE MAIN GEOMAGNETIC FIELD (IGRF MODEL)
% IF THE STANDARD GSM AND/OR SM COORDINATES ARE INTENDED TO BE USED, THEN SET VGSEX=-400.0 AND VGSEY=VGSEZ=0. IN THIS CASE, THE GSW COORDINATE SYSTEM BECOMES IDENTICAL TO THE STANDARD GSM.
recalc = @(y,d,h,m,s,vx,vy,vz)igrfmex(1,y,d,h,m,s,vx,vy,vz);

% SUBROUTINE IGRF_GSW_08 (XGSW,YGSW,ZGSW,HXGSW,HYGSW,HZGSW)
% INPUT PARAMETERS: XGSW,YGSW,ZGSW - CARTESIAN GEOCENTRIC SOLAR-WIND COORDINATES (IN UNITS RE=6371.2 KM)
igrfgsw = @(x,y,z)igrfmex(2,x,y,z); 

smgsw = @(x,y,z,d)igrfmex(3,x,y,z,d);
magsm = @(x,y,z,d)igrfmex(4,x,y,z,d);
geigeo = @(x,y,z,d)igrfmex(5,x,y,z,d);
geomag = @(x,y,z,d)igrfmex(6,x,y,z,d);
gswgse = @(x,y,z,d)igrfmex(7,x,y,z,d);
geogsw = @(x,y,z,d)igrfmex(8,x,y,z,d);