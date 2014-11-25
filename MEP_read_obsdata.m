clear;clc
indir = '../ENV-MA11_project_obsdata';
outdir = '../ENV-MA11_project_pics';
%% Read data
%% Meteorological data
fid = fopen([indir,'/','uea-wao-met_weybourne_20040101_met.txt'],'r');
meteo = textscan(fid, '%f %f %d %f %d %f %d %f %d %f %d %f %d %f %d', 'headerLines', 33);
fclose(fid);
meteo{2}(meteo{2}==9999.999) = NaN;
meteo{2}(meteo{3}==3) = NaN;
meteo{4}(meteo{4}==9999.999) = NaN;
meteo{4}(meteo{5}==3) = NaN;
meteo{6}(meteo{6}==9999.999) = NaN;
meteo{6}(meteo{7}==3) = NaN;
meteo{8}(meteo{8}==9999.999) = NaN;
meteo{8}(meteo{9}==3) = NaN;
meteo{10}(meteo{10}==9999.999) = NaN;
meteo{10}(meteo{11}==3) = NaN;
meteo{12}(meteo{12}==9999.999) = NaN;
meteo{12}(meteo{13}==3) = NaN;
meteo{14}(meteo{14}==9999.999) = NaN;
meteo{14}(meteo{15}==3) = NaN;
hur = meteo{2}; % relative humidity
ta = meteo{4}; % air temperature [C]
irrad = meteo{6}; % irradiance
netirrad = meteo{8}; % net irradiance
wspeed = meteo{10}; % wind speed [m s-1]
wdir = meteo{12}; % wind direction [deg]
ps = meteo{14}; % atmospheric pressure [mbar]
%% UEA O3 49C Analyser
fid = fopen([indir,'/','uea-o3_weybourne_20040419_ver1.txt'],'r');
O3_raw = textscan(fid, '%f %f %d', 'headerLines', 21);
fclose(fid);
O3_raw{2}(O3_raw{2}==9999) = NaN;
O3 = O3_raw{2};
%% UEA CraNOx
fid = fopen([indir,'/','uea-cranox_weybourne_20040419_ver1.txt'],'r');
cranox = textscan(fid, '%f %f %f', 'headerLines', 21);
fclose(fid);
cranox{2}(cranox{2}==9999) = NaN;
cranox{3}(cranox{3}==9999) = NaN;
NO = cranox{2};
NO2 = cranox{3};
%% Filter Radiometer J(NO2)
fid = fopen([indir,'/','leic_fr_jno2_weybourne_20040421.txt'],'r');
jNO2_raw = textscan(fid, '%f %f %d', 'headerLines', 26);
fclose(fid);
jNO2_raw{2}(jNO2_raw{2}==9999) = NaN;
jNO2_raw{2}(jNO2_raw{3}~=0) = NaN;
jNO2 = jNO2_raw{2};
%% [year,month,day,hour,min,sec,dayweek]
[year,month,day,hour,min,sec,dayweek] = julian2greg(116.9548611);