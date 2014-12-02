function [ta,ps,O3,NO2,NO,jNO2] = MEP_get_obsdata(iday)
indir = '../ENV-MA11_project_obsdata';
fid = fopen([indir,'/','uea-wao-met_weybourne_20040101_met.txt'],'r');
meteo = textscan(fid, '%f %f %d %f %d %f %d %f %d %f %d %f %d %f %d', 'headerLines', 33);
fclose(fid);
meteo{1} = meteo{1} - 1.5; % getting real dates
meteo{2}(meteo{2}==9999.999) = NaN; meteo{2}(meteo{3}==3) = NaN;
meteo{4}(meteo{4}==9999.999) = NaN; meteo{4}(meteo{5}==3) = NaN;
meteo{6}(meteo{6}==9999.999) = NaN; meteo{6}(meteo{7}==3) = NaN;
meteo{8}(meteo{8}==9999.999) = NaN; meteo{8}(meteo{9}==3) = NaN;
meteo{10}(meteo{10}==9999.999) = NaN; meteo{10}(meteo{11}==3) = NaN;
meteo{12}(meteo{12}==9999.999) = NaN; meteo{12}(meteo{13}==3) = NaN;
meteo{14}(meteo{14}==9999.999) = NaN; meteo{14}(meteo{15}==3) = NaN;
%% Read UEA O3 49C Analyser data
fid = fopen([indir,'/','uea-o3_weybourne_20040419_ver1.txt'],'r');
O3_raw = textscan(fid, '%f %f %d', 'headerLines', 21);
fclose(fid);
O3_raw{1} = O3_raw{1} - 1.5; % getting real dates
O3_raw{2}(O3_raw{2}==9999) = NaN;
%% Read UEA CraNOx data
fid = fopen([indir,'/','uea-cranox_weybourne_20040419_ver1.txt'],'r');
cranox = textscan(fid, '%f %f %f', 'headerLines', 21);
fclose(fid);
cranox{1} = cranox{1} - 1.5; % getting real dates
cranox{2}(cranox{2}==9999 | cranox{2} < 0) = NaN;
cranox{3}(cranox{3}==9999 | cranox{3} < 0) = NaN;
%% Read Filter Radiometer j(NO2) data
fid = fopen([indir,'/','leic_fr_jno2_weybourne_20040421.txt'],'r');
jNO2_raw = textscan(fid, '%f %f %d', 'headerLines', 26);
fclose(fid);
jNO2_raw{1} = jNO2_raw{1} - 1.5; % getting real dates
jNO2_raw{2}(jNO2_raw{2}==9999) = NaN;
jNO2_raw{2}(jNO2_raw{3}~=0) = NaN;
%% Select day
ndays = 38; % total number of days of observations
% Set dates array
dates(1,1:9) = 22:30; dates(2,1:9) = 4;
dates(1,10:38) = 1:29; dates(2,10:38) = 5;
% Start dates indices
ibegc1 = 1441; % for jNO2
ibegc2 = 5108; % for NO2, NO, O3
ibegm = 16129; % 16129 22 Apr 00:00:00 (111.5)

hur = meteo{2}((ibegm+144*(iday-1)):(ibegm+143*iday)); % relative humidity
ta = meteo{4}((ibegm+144*(iday-1)):(ibegm+143*iday)); % air temperature [C]
irrad = meteo{6}((ibegm+144*(iday-1)):(ibegm+143*iday)); % irradiance
netirrad = meteo{8}((ibegm+144*(iday-1)):(ibegm+143*iday)); % net irradiance
wspeed = meteo{10}((ibegm+144*(iday-1)):(ibegm+143*iday)); % wind speed [m s-1]
wdir = meteo{12}((ibegm+144*(iday-1)):(ibegm+143*iday)); % wind direction [deg]
ps = meteo{14}((ibegm+144*(iday-1)):(ibegm+143*iday)); % atmospheric pressure [mbar]

NO = cranox{2}((ibegc2+1440*(iday-1)):(ibegc2+1439*iday));
NO2 = cranox{3}((ibegc2+1440*(iday-1)):(ibegc2+1439*iday));
jNO2 = jNO2_raw{2}((ibegc1+1440*(iday-1)):(ibegc1+1439*iday));
O3 = O3_raw{2}((ibegc2+1440*(iday-1)):(ibegc2+1439*iday));
end