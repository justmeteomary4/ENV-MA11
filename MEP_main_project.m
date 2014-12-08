clear; close all; clc
outdir = '../ENV-MA11_project_pics';
%% Model parameters
kB = 1.380658*1e-19; % Boltzmann's constant [cm-3 hPa K-1 molec-1]
hbeg = 3;
hend = 20;
hstep = 0.25; % 0.25 - 15 min, 0.50 - 30 min, 1 - 1 hour
hsteps = hbeg:hstep:hend;
it = 0; % t0 loop
for tstart = hsteps
it = it + 1
t0 = 60*60*tstart; % time of the day [s]
nt = 60*60*0.25; % simulation time [s]
dt = 10; % time step [s]
nsteps = numel(t0:dt:nt); % number of time steps
SCHEME = 2; % 1 - Forward Euler (nt = 1 s, dt = 1e-7 s); 2 - Backward Euler; 3 - MIE (nt = 2000 s, dt = 10 s)
%% Initial condtions
EXPERIMENT  = 2;
switch EXPERIMENT
    case 1 % test
        T = 298.; % air temperature [Kelvin], ~= 25 degrees Celsius
%         p = 1013.25; % air pressure [hPa], == 1 atm (standard atmosphere)
%         M = p/(kB*T);
%         O2 = 0.2095*M; % assuming O2 mixiting ratio = 0.2095 and the air is dry
        % Number densities (concentrations) [molec cm-3]
        O3 = 0.;
        O = 0.;
        NO2 = 1e10; % ~= 400 pptv (mixing ratio)
        NO = 1e12; % ~= 40 ppbv (mixing ratio)
        % Rate coefficients
        k1 = 1.4*1e3*exp(1175/T); % O -> O3, this coefficient has already been multiplied by O2 and M
        jNO2 = 1.7*1e-2; % NO2+hv -> NO+O, <420 nm
        k3 = 1.8*1e-12*exp(-1370/T); % NO+O3 -> NO2+O2
    case 2 % date
        date = '27Apr';
        [T,ps,O3_obs,NO2_obs,NO_obs,jNO2_obs] = MEP_get_obsdata(38);
        l = find(~isnan(jNO2_obs),1,'first'); % find first point with photolysis (121st), jNO2 counter
        l = max(ceil(t0/60),l); % chose bigger l index with photolysis
        nt = min(t0+nt,(find(~isnan(jNO2_obs),1,'last'))*60); % the last existing point with photolysis 82740, (61 min from the end)
        while mod(l,10) ~= 0
            l = l + 1;
        end
        k = l/10; % temperature counter
        l = l - 1;
        k = k - 1;
        t0 = max(t0,l*60);
        % Mixing ratios [ppbv]
        O3 = O3_obs(l); %disp(['O3 = ' num2str(O3)])
        O = 0.; %disp(['O = ' num2str(O)])
        NO2 = NO2_obs(l); %disp(['NO2 = ' num2str(NO2)])
        NO = NO_obs(l); %disp(['NO = ' num2str(NO)])
end
%% Numerical solution
switch SCHEME
    case 1 % Forward Euler
        % Preallocating space for variables
        O3forw(1:nsteps) = NaN; Oforw(1:nsteps) = NaN; NOforw(1:nsteps) = NaN; NO2forw(1:nsteps) = NaN;
        O3forw(1) = O3(1); Oforw(1) = O(1); NOforw(1) = NO(1); NO2forw(1) = NO2(1);
        i = 0;
        for t = t0:dt:nt
            i = i + 1;
            NO2forw(i+1) = NO2forw(i)+dt*(k3*NOforw(i)*O3forw(i) - jNO2*NO2forw(i));
            NOforw(i+1) = NOforw(i)+dt*(jNO2*NO2forw(i) - k3*NOforw(i)*O3forw(i));
            Oforw(i+1) = Oforw(i) + dt*(jNO2*NO2forw(i) - k1*Oforw(i));
            O3forw(i+1) = O3forw(i)+dt*(k1*Oforw(i) - k3*NOforw(i)*O3forw(i));
        end
    case 2 % Backward Euler
            % Preallocating space for variables
            O3back(1:nsteps) = NaN; Oback(1:nsteps) = NaN; NOback(1:nsteps) = NaN; NO2back(1:nsteps) = NaN;
            O3back(1) = O3; Oback(1) = O; NOback(1) = NO; NO2back(1) = NO2;
            switch EXPERIMENT
                case 1 % test
                    i = 0;
                    for t = t0:dt:nt
                        i = i+1;
                        NO2back(i+1) = (NO2back(i) + dt*k3*NOback(i)*O3back(i))/(1 + dt*jNO2);
                        NOback(i+1) = (NOback(i) + dt*jNO2*NO2back(i))/(1 + dt*k3*O3back(i));
                        Oback(i+1) = (Oback(i) + dt*jNO2*NO2back(i))/(1 + dt*k1);
                        O3back(i+1) = (O3back(i) + dt*k1*Oback(i))/(1 + dt*k3*NOback(i));
                    end
                case 2 % 27 May
                    i = 0; % iteration
                    atimesec = (t0+dt):dt:nt; % absolute time of the day [s]
                    for t = atimesec
                        if mod(i*10,600)==0 % update temperature every 10 min
                            k = k+1;
                            k1 = 1.4*1e3*exp(1175/T(k));
                            k3 = 1.8*1e-12*exp(-1370/T(k));
                        end
                        if mod(i*10,60)==0 % update jNO2 every 1 min
                            l = l+1;
                            jNO2 = jNO2_obs(l);
                        end
                        i = i+1;
                        NO2back(i+1) = (NO2back(i) + dt*k3*NOback(i)*O3back(i))/(1 + dt*jNO2);
                        NOback(i+1) = (NOback(i) + dt*jNO2*NO2back(i))/(1 + dt*k3*O3back(i));
                        Oback(i+1) = (Oback(i) + dt*jNO2*NO2back(i))/(1 + dt*k1);
                        O3back(i+1) = (O3back(i) + dt*k1*Oback(i))/(1 + dt*k3*NOback(i));
                    end
                    tmodel(it) = nt;
                    NO2backlast(it) = NO2back(i+1);
                    NObacklast(it) = NOback(i+1);
                    Obacklast(it) = Oback(i+1);
                    O3backlast(it) = O3back(i+1);                    
            end
    case 3 % MIE
        Np = 30; % limit of number of iterations; 5 for large set of reastions, 30-50 for small set of reactions
        itermax = 1000; % working limit of number of iterations
        % Preallocating space for variables
        Nmie(1:nsteps,1:4) = NaN;
        P(1:nsteps,1:4) = NaN;
        L(1:nsteps,1:4) = NaN;
        % Step 1
        % Initialize Backward Euler concentrations and maximum Backward Euler concentrations for each
        % active species with concentrations from the beginning of the simulation
        Nmiemax = [NO2 NO O O3];
        Nmieback = [NO2 NO O O3];
        Nmie(1,:) =  [NO2 NO O O3];
        % Step 3
        % Estimate production rates and implicit loss coefficients for each species
        % Time loop
        for i = 2:nsteps
%             disp(['i = ' num2str(i) '; Time: ' num2str(t0+dt*(i-1)) 's'])
            % Iteration loop
            m = 0;
            np = 0;
            while np <= Np
                m = m + 1;
                if m > itermax
                    error(['Iterations do not converge for the time step ' num2str(i)])
                end
                % Production rates for each species
                P(i,1) = k3*Nmieback(2)*Nmieback(4); % NO2
                P(i,2) = jNO2*Nmieback(1); % NO
                P(i,3) = jNO2*Nmieback(1); % O
                P(i,4) = k1*Nmieback(3); % O3
                % Implicit loss coef?cients for each species
                L(i,1) = jNO2; % NO2
                L(i,2) = k3*Nmieback(4); % NO
                L(i,3) = k1; % O
                L(i,4) = k3*Nmieback(2); % O3
                % Step 4 and 5
                % Calculate backward and forward Euler concentrations for all species at iteration m +1
                for s = 1:4
                    Nmieback(s) = (Nmie(i-1,s) + dt*P(i,s))/(1 + dt*L(i,s));
                    Nmieforw(s) = Nmie(i-1,s) + dt*(P(i,s) - L(i,s));
                end
                % Step 6
                % Check convergence
                if all(Nmieforw(:) >= 0)
                    np =np + 1;
                else
                    np = 0;
                end
                work = Nmieback; % working variable
                for s = 1:4
                    % Limit current backward Euler concentrations
                    Nmieback(s) = min(Nmieback(s), Nmiemax(s));
                    % Update maximum backward Euler concentrations
                    Nmiemax(s) = max(work(s), Nmie(i-1,s));
                end
            end
            % Step 7
            % The final concentrations
            Nmie(i,:) = Nmieback(:);
        end
end
end % t0 loop
%% Plot results
switch SCHEME
    case 1 % Forward Euler
        switch EXPERIMENT
            case 1 % test
                figure; xend = 1e7;
                subplot(2,2,1); plot(NO2forw(1:end-1),'ko','MarkerSize',3); title('NO2');
                xlabel('seconds'); ylabel('molecules cm-3'); set(gca,'Xtick',linspace(0,1e7,3),'XTickLabel',0:0.5:1); xlim([0 xend])
                subplot(2,2,2); plot(NOforw(1:end-1),'ks','MarkerSize',3); title('NO');
                xlabel('seconds'); ylabel('molecules cm-3'); set(gca,'Xtick',linspace(0,1e7,3),'XTickLabel',0:0.5:1); xlim([0 xend])
                subplot(2,2,3); plot(Oforw(1:end-1),'k^','MarkerSize',3); title('O');
                xlabel('seconds'); ylabel('molecules cm-3'); set(gca,'Xtick',linspace(0,1e7,3),'XTickLabel',0:0.5:1); xlim([0 xend])
                subplot(2,2,4); plot(O3forw(1:end-1),'kd','MarkerSize',3); title('O3');
                xlabel('seconds'); ylabel('molecules cm-3'); set(gca,'Xtick',linspace(0,1e7,3),'XTickLabel',0:0.5:1); xlim([0 xend])
                imgname= strcat(outdir,'/pics_test/','test_forwEuler','.png');
%                 set(gcf,'visible','off')
%                 print(gcf,'-dpng','-r300',imgname);
        end
    case 2 % Backward Euler
        switch EXPERIMENT
            case 1 % test
                figure;
                subplot(2,2,1); plot(NO2back(1:end-1),'bd','MarkerSize',3); title('NO2');
                xlabel('seconds'); ylabel('molecules cm-3'); set(gca,'Xtick',0:50:200,'XTickLabel',0:500:2000); xlim([0 nsteps-1]);
                subplot(2,2,2); plot(NOback(1:end-1),'b^','MarkerSize',3); title('NO');
                xlabel('seconds'); ylabel('molecules cm-3'); set(gca,'Xtick',0:50:200,'XTickLabel',0:500:2000); xlim([0 nsteps-1]);
                subplot(2,2,3); plot(Oback(1:end-1),'bs','MarkerSize',3); title('O');
                xlabel('seconds'); ylabel('molecules cm-3'); set(gca,'Xtick',0:50:200,'XTickLabel',0:500:2000); xlim([0 nsteps-1]);
                subplot(2,2,4); plot(O3back(1:end-1),'bo','MarkerSize',3); title('O3');
                xlabel('seconds'); ylabel('molecules cm-3'); set(gca,'Xtick',0:50:200,'XTickLabel',0:500:2000); xlim([0 nsteps-1]);
                imgname= strcat(outdir,'/pics_test/','test_backEuler','.png');
%                 set(gcf,'visible','off')
%                 print(gcf,'-dpng','-r300',imgname);
            case 2 % 27 May
                figure(1); % plot model data
                subplot(2,2,1); plot(NO2back(1:end-1),'bd','MarkerSize',3); title('NO2');
                xlabel('seconds'); ylabel('ppbv'); xlim([1 numel(atimesec)]);
                set(gca,'XTick',[1 numel(atimesec)],'Xticklabel',[atimesec(1) atimesec(end)]); 
                subplot(2,2,2); plot(NOback(1:end-1),'b^','MarkerSize',3); title('NO');
                xlabel('seconds'); ylabel('ppbv'); xlim([1 numel(atimesec)]);
                set(gca,'XTick',[1 numel(atimesec)],'Xticklabel',[atimesec(1) atimesec(end)]);
                subplot(2,2,3); plot(Oback(1:end-1),'bs','MarkerSize',3); title('O');
                xlabel('seconds'); ylabel('ppbv'); xlim([1 numel(atimesec)]);
                set(gca,'XTick',[1 numel(atimesec)],'Xticklabel',[atimesec(1) atimesec(end)]);
                subplot(2,2,4); plot(O3back(1:end-1),'bo','MarkerSize',3); title('O3');
                xlabel('seconds'); ylabel('ppbv'); xlim([1 numel(atimesec)]);
                set(gca,'XTick',[1 numel(atimesec)],'Xticklabel',[atimesec(1) atimesec(end)]);
                imgname= strcat(outdir,'/pics_model&obs/',date,'_backEuler_',num2str(hstep),'_part','.png');
%                 set(gcf,'visible','off')
%                 print(gcf,'-dpng','-r300',imgname);

                figure(2); % overlay model and observational data
                subplot(2,2,1); plot(tmodel,NO2backlast,'bd','MarkerSize',3); title('NO2'); hold on;
                plot(hsteps*60*60,NO2_obs(hbeg*60:hstep*60:hend*60),'r-')
                xlabel('seconds'); ylabel('ppbv');
                ylim([min(min(NO2backlast(:)),min(NO2_obs(:))) max(max(NO2backlast(:),max(NO2_obs(:))))])
                subplot(2,2,2); plot(tmodel,NObacklast,'b^','MarkerSize',3); title('NO'); hold on;
                plot(hsteps*60*60,NO_obs(hbeg*60:hstep*60:hend*60),'r-')
                xlabel('seconds'); ylabel('ppbv');
                ylim([min(min(NObacklast(:)),min(NO_obs(:))) max(max(NObacklast(:),max(NO_obs(:))))])
                subplot(2,2,3); plot(tmodel,Obacklast,'bs','MarkerSize',3); title('O'); 
                xlabel('seconds'); ylabel('ppbv');
                ylim([min(Obacklast(:)) max(Obacklast(:))])
                subplot(2,2,4); plot(tmodel,O3backlast,'bo','MarkerSize',3); title('O3'); hold on;
                plot(hsteps*60*60,O3_obs(hbeg*60:hstep*60:hend*60),'r-')
                xlabel('seconds'); ylabel('ppbv');
                ylim([min(min(O3backlast(:)),min(O3_obs(:))) max(max(O3backlast(:),max(O3_obs(:))))])
                imgname= strcat(outdir,'/pics_model&obs/',date,'_backEuler_vs_obs_',num2str(hstep),'_whole','.png');
%                 set(gcf,'visible','off')
                print(gcf,'-dpng','-r300',imgname);
        end
    case 3 % MIE
        switch EXPERIMENT
            case 1 %test
                figure;
                subplot(2,2,1); plot(Nmie(1:end-1,1),'rd','MarkerSize',3); title('NO2'); hold on; plot(NO2back(1:end-1),'bd','MarkerSize',3);
                xlabel('seconds'); ylabel('molecules cm-3'); set(gca,'Xtick',0:50:200,'XTickLabel',0:500:2000); xlim([0 nsteps-1]);
                subplot(2,2,2); plot(Nmie(1:end-1,2),'r^','MarkerSize',3); title('NO'); hold on; plot(NOback(1:end-1),'b^','MarkerSize',3);
                xlabel('seconds'); ylabel('molecules cm-3'); set(gca,'Xtick',0:50:200,'XTickLabel',0:500:2000); xlim([0 nsteps-1]);
                subplot(2,2,3); plot(Nmie(1:end-1,3),'rs','MarkerSize',3); title('O'); hold on; plot(Oback(1:end-1),'bs','MarkerSize',3);
                xlabel('seconds'); ylabel('molecules cm-3'); set(gca,'Xtick',0:50:200,'XTickLabel',0:500:2000); xlim([0 nsteps-1]);
                subplot(2,2,4); plot(Nmie(1:end-1,4),'ro','MarkerSize',3); title('O3'); hold on; plot(O3back(1:end-1),'bo','MarkerSize',3);
                xlabel('seconds'); ylabel('molecules cm-3'); set(gca,'Xtick',0:50:200,'XTickLabel',0:500:2000); xlim([0 nsteps-1]);
                imgname= strcat(outdir,'/pics_test/','test_mie&backEuler','.png');
%                 set(gcf,'visible','off')
%                 print(gcf,'-dpng','-r300',imgname);
        end
end