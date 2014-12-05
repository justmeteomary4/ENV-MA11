clear; close all; clc
%% Set input and output directories
indir = '../ENV-MA11_project_obsdata';
outdir = '../ENV-MA11_project_pics';
%% Read observational data
% In source files time (in julian days) is greater than real time by 1.5, year starts from 1 rather than from -0.5. 
%% Read meteorological data
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
hur = meteo{2}; % relative humidity
ta = meteo{4}; % air temperature [C]
irrad = meteo{6}; % irradiance
netirrad = meteo{8}; % net irradiance
wspeed = meteo{10}; % wind speed [m s-1]
wdir = meteo{12}; % wind direction [deg]
ps = meteo{14}; % atmospheric pressure [mbar]
%% Read UEA O3 49C Analyser data
fid = fopen([indir,'/','uea-o3_weybourne_20040419_ver1.txt'],'r');
O3_raw = textscan(fid, '%f %f %d', 'headerLines', 21);
fclose(fid);
O3_raw{1} = O3_raw{1} - 1.5; % getting real dates
O3_raw{2}(O3_raw{2}==9999) = NaN;
O3 = O3_raw{2};
%% Read UEA CraNOx data
fid = fopen([indir,'/','uea-cranox_weybourne_20040419_ver1.txt'],'r');
cranox = textscan(fid, '%f %f %f', 'headerLines', 21);
fclose(fid);
cranox{1} = cranox{1} - 1.5; % getting real dates
cranox{2}(cranox{2}==9999 | cranox{2} < 0) = NaN;
cranox{3}(cranox{3}==9999 | cranox{3} < 0) = NaN;
NO = cranox{2};
NO2 = cranox{3};
%% Read Filter Radiometer j(NO2) data
fid = fopen([indir,'/','leic_fr_jno2_weybourne_20040421.txt'],'r');
jNO2_raw = textscan(fid, '%f %f %d', 'headerLines', 26);
fclose(fid);
jNO2_raw{1} = jNO2_raw{1} - 1.5; % getting real dates
jNO2_raw{2}(jNO2_raw{2}==9999) = NaN;
jNO2_raw{2}(jNO2_raw{3}~=0) = NaN;
jNO2 = jNO2_raw{2};
%% julian2greg convector
%[year,month,day,hour,min,sec,dayweek] = julian2greg(116.9548611);
%% Plot observational data
% Time of the day in source files: night - 0.0, day - 0.5. Time of the day here (real): night - 0.5, day - 0.0
% Start/end dates and their indices:
% Variable: start date time (real julian day) [index] / end date time (real julian day) [index]
% jNO2: 22 Apr 00:00:00 (111.5) [1441] / 29 May 23:59:01 (149.4993056) [56160]
% NO2, NO, O3: 21 Apr 23:59:59 (111.499988) [5108] / 29 May 23:59:01 (149.499306) [59826]
ndays = 38; % total number of days of observations
% Set dates array
dates(1,1:9) = 22:30; dates(2,1:9) = 4;
dates(1,10:38) = 1:29; dates(2,10:38) = 5;
% Start dates indices
ibegc1 = 1441; % for jNO2
ibegc2 = 5108; % for NO2, NO, O3
ibegm = 16129; % 16129 22 Apr 00:00:00 (111.5)
cvec = jet(ndays); % colormap
plottype = 2;
switch plottype
    case 1
        figure(1); % Chemistry for the whole period on one graph
        subplot(2,2,1); plot(NO2(5108:59826)); title('NO2'); 
        subplot(2,2,2); plot(NO(5108:59826)); title('NO');
        subplot(2,2,3); plot(jNO2(1441:56160)); title('jNO2');
        subplot(2,2,4); plot(O3(5108:59826)); title('O3');
    case 2
        for i = 1:ndays
            figure(109);
            subplot(2,2,1); plot(NO2((ibegc2+1440*(i-1)):(ibegc2+1439*i)),'color',cvec(i,:)); %hold on;
            title('NO2'); xlabel('min'); ylabel('ppbv'); xlim([0 1440]); ylim([0 max(NO2(:))]);
            subplot(2,2,2); plot(NO((ibegc2+1440*(i-1)):(ibegc2+1439*i)),'color',cvec(i,:)); %hold on;
            title('NO'); xlabel('min'); ylabel('ppbv'); xlim([0 1440]); ylim([0 max(NO(:))]);
            subplot(2,2,3); plot(jNO2((ibegc1+1440*(i-1)):(ibegc1+1439*i)),'color',cvec(i,:)); %hold on;
            title('jNO2'); xlabel('min'); ylabel('1/s'); xlim([0 1440]); ylim([0 max(jNO2(:))]); 
            subplot(2,2,4); plot(O3((ibegc2+1440*(i-1)):(ibegc2+1439*i)),'color',cvec(i,:)); %hold on;
            title('O3'); xlabel('min'); ylabel('ppbv'); xlim([0 1440]); ylim([0 max(O3(:))]);
            imgname= strcat(outdir,'/pics_obsdata_chem/','obsdata_chem_',datestr([2004,dates(2,i),dates(1,i),0,0,0],'mmm dd'),'.png');
            suptitle(datestr([2004,dates(2,i),dates(1,i),0,0,0],'mmm dd'));
            set(gcf,'visible','off')
            print(gcf,'-dpng','-r300',imgname);
            
            figure(110);
            subplot(2,2,1); plot(ta((ibegm+144*(i-1)):(ibegm+143*i)),'color',cvec(i,:));
            title('Air temperature'); xlabel('min'); ylabel('degrees Celsius'); xlim([0 144]); ylim([min(ta(:)) max(ta(:))]);
            subplot(2,2,2); plot(ps((ibegm+144*(i-1)):(ibegm+143*i)),'color',cvec(i,:));
            title('Surface air pressure'); xlabel('min'); ylabel('hPa'); xlim([0 144]); ylim([min(ps(:)) max(ps(:))]);
            subplot(2,2,3); plot(wdir((ibegm+144*(i-1)):(ibegm+143*i)),'color',cvec(i,:));
            title('Wind direction'); xlabel('min'); ylabel('degrees'); xlim([0 144]); ylim([0 360]);
            subplot(2,2,4); plot(wspeed((ibegm+144*(i-1)):(ibegm+143*i)),'color',cvec(i,:));
            title('Wind speed'); xlabel('min'); ylabel('m/s'); xlim([0 144]); ylim([min(wspeed(:)) max(wspeed(:))]);
            imgname= strcat(outdir,'/pics_obsdata_meteo/','obsdata_meteo_',datestr([2004,dates(2,i),dates(1,i),0,0,0],'mmm dd'),'.png');
            suptitle(datestr([2004,dates(2,i),dates(1,i),0,0,0],'mmm dd'));
            set(gcf,'visible','off')
            print(gcf,'-dpng','-r300',imgname);
        end
end