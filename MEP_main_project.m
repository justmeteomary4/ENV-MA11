%clear;clc;close
%% Model parameters
T = 298.; % absolute temperature [Kelvins], +25 C
p = 1013.; % total air pressure [hPa]
kB = 1.380658*1e-19; % Boltzmann's constant [cm-3 hPa K-1 molec-1]
t0 = 0;
nt = 2000; % simulation time [s]
dt = 10; % time step [s]
nsteps = numel(t0:dt:nt); % number of time steps
scheme = 3; % 1 - Forward Euler (nt = 1 s, dt = 1e-4) s); 2 - Backward Euler (nt = 2000 s, dt = 10 s); 3 - MIE
Np = 30; % 5 for large set of reastions, 30-50 for small set of reactions
itermax = 100; % empirically derived reosanoble ilimit of number of iterations
%% Initial concentrations
O3 = 0.; % molec cm-3
M = p/(kB*T);
O2 = 0.2095*M; % assuming mixiting ratio of O2 = 0.2095 and the air is dry
O = 0.;
NO = 1e12;
NO2 = 1e10;
%% Rate coefficients
k1 = 1.4*1e3*exp(1175/T); % O -> O3, this coefficient has already been multiplied by O2 and M
J = 1.7*1e-2; % NO2+hv -> NO+O, <420 nm
k3 = 1.8*1e-12*exp(-1370/T); % NO+O3 -> NO2+O2
%% Numerical solution
switch scheme
    case 1 % Forward Euler
        % Preallocating space for variables
        O3forw(1:nsteps) = NaN; Oforw(1:nsteps) = NaN; NOforw(1:nsteps) = NaN; NO2forw(1:nsteps) = NaN;
        O3forw(1) = O3(1); Oforw(1) = O(1); NOforw(1) = NO(1); NO2forw(1) = NO2(1);
        i = 0;
%         for t = t0:dt:nt % taking O2 and M concentrations into account
%             i = i + 1;
%             NO2forw(i+1) = NO2forw(i)+dt*(k1*NOforw(i)*O3forw(i) - J*NO2forw(i));
%             NOforw(i+1) = NOforw(i)+dt*(J*NO2forw(i) - k1*NOforw(i)*O3forw(i));
%             Oforw(i+1) = Oforw(i) + dt*(J*NO2forw(i) - k3*Oforw(i)*O2*M);
%             O3forw(i+1) = O3forw(i)+dt*(k3*Oforw(i)*O2*M - k1*NOforw(i)*O3forw(i));
%         end
        for t = t0:dt:nt % do not taking O2 and M concentrations into account
            i = i + 1;
            NO2forw(i+1) = NO2forw(i)+dt*(k3*NOforw(i)*O3forw(i) - J*NO2forw(i));
            NOforw(i+1) = NOforw(i)+dt*(J*NO2forw(i) - k3*NOforw(i)*O3forw(i));
            Oforw(i+1) = Oforw(i) + dt*(J*NO2forw(i) - k1*Oforw(i));
            O3forw(i+1) = O3forw(i)+dt*(k1*Oforw(i) - k3*NOforw(i)*O3forw(i));
        end
    case 2 % Backward Euler
            % Preallocating space for variables
            O3back(1:nsteps) = NaN; Oback(1:nsteps) = NaN; NOback(1:nsteps) = NaN; NO2back(1:nsteps) = NaN;
            O3back(1) = O3(1); Oback(1) = O(1); NOback(1) = NO(1); NO2back(1) = NO2(1);
            i = 0;
            for t = t0:dt:nt
                i = i+1;
                NO2back(i+1) = (NO2back(i) + dt*k3*NOback(i)*O3back(i))/(1 + dt*J);
                NOback(i+1) = (NOback(i) + dt*J*NO2back(i))/(1 + dt*k3*O3back(i));
                Oback(i+1) = (Oback(i) + dt*J*NO2back(i))/(1 + dt*k1);
                O3back(i+1) = (O3back(i) + dt*k1*Oback(i))/(1 + dt*k3*NOback(i));
            end
    case 3 % Multistep implicit–explicit (MIE)
        % Preallocating space for variables
        Nmie(1:nsteps,1:4) = NaN;
        P(1:nsteps,1:4) = NaN;
        L(1:nsteps,1:4) = NaN;
        %% Step 1
        % Initialize Backward Euler concentrations and maximum Backward Euler concentrations for each
        % active species with concentrations from the beginning of the simulation
        Nmiemax = [NO2 NO O O3];
        Nmieback = [NO2 NO O O3];
        Nmie(1,:) =  [NO2 NO O O3];
        %% Step 2
        % Estimate reaction rates by multiplying rate coeffcients by Backward Euler concentrations
        % Skipped, go directly to production and loss rates
%         R1 = k1*NOmieback(1)*O3mieback(1);
%         R2 = J*NO2mieback(1);
%         R3 = k3*Omieback(1);
        %% Step 3
        % Estimate production rates and implicit loss coef?cients for each species
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
                P(i,2) = J*Nmieback(1); % NO
                P(i,3) = J*Nmieback(1); % O
                P(i,4) = k1*Nmieback(3); % O3
                % Implicit loss coef?cients for each species
                L(i,1) = J; % NO2
                L(i,2) = k3*Nmieback(4); % NO
                L(i,3) = k1; % O
                L(i,4) = k3*Nmieback(2); % O3
                %% Step 4 and 5
                % Calculate backward and forward Euler concentrations for all species at iteration m +1
                for s = 1:4
                    Nmieback(s) = (Nmie(i-1,s) + dt*P(i,s))/(1 + dt*L(i,s));
                    Nmieforw(s) = Nmie(i-1,s) + dt*(P(i,s) - L(i,s));
                end
                %% Step 6
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
            %% Step 7
            % The final concentrations
            Nmie(i,:) = Nmieback(:);
        end
end
%% Plot results
switch scheme
    case 1
        figure;
        end1 = 20;
        subplot(2,2,1); plot(O3forw(1:end1),'k'); title('O3')
        subplot(2,2,2); plot(Oforw(1:end1),'k'); title('O')
        subplot(2,2,3); plot(NOforw(1:end1),'k'); title('NO')
        subplot(2,2,4); plot(NO2forw(1:end1),'k'); title('NO2')
        %subplot(2,2,1); plot(t0:dt:nt,O3forw(1:end-1),'b'); title('O3')
        %subplot(2,2,2); plot(t0:dt:nt,Oforw(1:end-1),'b'); title('O')
        %subplot(2,2,3); plot(t0:dt:nt,NOforw(1:end-1),'b'); title('NO')
        %subplot(2,2,4); plot(t0:dt:nt,NO2forw(1:end-1),'b'); title('NO2')
    case 2
        figure;
        subplot(2,2,1); plot(NO2back(1:end-1),'b'); title('NO2')
        subplot(2,2,2); plot(NOback(1:end-1),'b'); title('NO')
        subplot(2,2,3); plot(Oback(1:end-1),'b'); title('O')
        subplot(2,2,4); plot(O3back(1:end-1),'b'); title('O3')
    case 3
        figure;
        subplot(2,2,1); plot(Nmie(1:end-1,1),'r'); title('NO2'); hold on; plot(NO2back(1:end-1),'b'); title('NO2')
        subplot(2,2,2); plot(Nmie(1:end-1,2),'r'); title('NO'); hold on; plot(NOback(1:end-1),'b'); title('NO')
        subplot(2,2,3); plot(Nmie(1:end-1,3),'r'); title('O'); hold on; plot(Oback(1:end-1),'b'); title('O')
        subplot(2,2,4); plot(Nmie(1:end-1,4),'r'); title('O3'); hold on; plot(O3back(1:end-1),'b'); title('O3')
end