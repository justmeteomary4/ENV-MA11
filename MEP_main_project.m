clear;clc;close
%% Model parameters
t0 = 0;
nt = 400; % simulation time
dt = 10; % time step [s]
nsteps = numel(t0:dt:nt);
T = 298; % absolute temperature [Kelvins], +25 C
p = 920; % total air pressure [hPa]
kB = 1.380658*10^(-19); % Boltzmann's constant [cm-3 hPa K-1 molec-1]
scheme = 3; % 1 - Forward Euler (nt = 1 s, dt = 10^(-4) s); 2 - Backward Euler (nt = 2000 s, dt = 10 s); 3 - MIE
Np = 30; % 5 for large set of reastions, 30-50 for small set of reactions
%% Initial concentrations
O3 = 0; % molec cm-3
M = p/(kB*T);
O2 = 0.2095*M; % assuming mixiting ratio of O2 = 0.2095 and air is dry
O = 0;
NO = 10^(12);
NO2 = 10^(10);
%% Rate coefficients
k1 = 1.8*10^(-12)*exp(-1370/T); % NO+O3 -> NO2+O2
J = 1.7*10^(-2);                        % NO2+hv -> NO+O, <420 nm
k3 = 1.4*10^(3)*exp(1175/T);     % O+O2+M -> O3+M; simply O -> O3
%% Numerical solution
switch scheme
    case 1 % Forward Euler
        % Preallocating space for variables
        O3forw(1:nsteps) = NaN; Oforw(1:nsteps) = NaN; NOforw(1:nsteps) = NaN; NO2forw(1:nsteps) = NaN;
        O3forw(1) = O3(1); Oforw(1) = O(1); NOforw(1) = NO(1); NO2forw(1) = NO2(1);
        i = 0;
%         for t = t0:dt:nt % taking O2 and M concentrations into account
%             i = i + 1;
%             NOforw(i+1) = NOforw(i)+dt*(J*NO2forw(i) - k1*NOforw(i)*O3forw(i));
%             NO2forw(i+1) = NO2forw(i)+dt*(k1*NOforw(i)*O3forw(i) - J*NO2forw(i));
%             Oforw(i+1) = Oforw(i) + dt*(J*NO2forw(i) - k3*Oforw(i)*O2*M);
%             O3forw(i+1) = O3forw(i)+dt*(k3*Oforw(i)*O2*M - k1*NOforw(i)*O3forw(i));
%         end
        for t = t0:dt:nt % do not taking O2 and M concentrations into account
            i = i + 1;
            NOforw(i+1) = NOforw(i)+dt*(J*NO2forw(i) - k1*NOforw(i)*O3forw(i));
            NO2forw(i+1) = NO2forw(i)+dt*(k1*NOforw(i)*O3forw(i) - J*NO2forw(i));
            Oforw(i+1) = Oforw(i) + dt*(J*NO2forw(i) - k3*Oforw(i));
            O3forw(i+1) = O3forw(i)+dt*(k3*Oforw(i) - k1*NOforw(i)*O3forw(i));
        end
        case 2 % Backward Euler
            % Preallocating space for variables
            O3back(1:nsteps) = NaN; Oback(1:nsteps) = NaN; NOback(1:nsteps) = NaN; NO2back(1:nsteps) = NaN;
            O3back(1) = O3(1); Oback(1) = O(1); NOback(1) = NO(1); NO2back(1) = NO2(1);
            i = 0;
            for t = t0:dt:nt
                i = i+1;
                NOback(i+1) = (NOback(i) + dt*J*NO2back(i))/(1 + dt*k1*O3back(i));
                NO2back(i+1) = (NO2back(i) + dt*k1*NOback(i)*O3back(i))/(1 + dt*J);
                Oback(i+1) = (Oback(i) + dt*J*NO2back(i))/(1 + dt*k3);
                O3back(i+1) = (O3back(i) + dt*k3*Oback(i))/(1 + dt*k1*NOback(i));
            end
    case 3 % MIE
        % Preallocating space for variables
        O3mie(1:nsteps) = NaN; Omie(1:nsteps) = NaN; NOmie(1:nsteps) = NaN; NO2mie(1:nsteps) = NaN;
%         O3mieback(1:nsteps) = NaN; Omieback(1:nsteps) = NaN; NOmieback(1:nsteps) = NaN; NO2mieback(1:nsteps) = NaN;
%         O3mieforw(1:nsteps) = NaN; Omieforw(1:nsteps) = NaN; NOmieforw(1:nsteps) = NaN; NO2mieforw(1:nsteps) = NaN;
        P(1:nsteps,1:4) = NaN;
        L(1:nsteps,1:4) = NaN;
        %% First stage
        % Initialize Backward Euler concentrations and maximum Backward Euler concentrations for each
        % active species with concentrations from the beginning of the simulation
        O3mie(1) = O3; Omie(1) = O; NOmie(1) = NO; NO2mie(1) = NO2;
        O3miemax = O3; Omiemax = O; NOmiemax = NO; NO2miemax = NO2;
        O3mieback = O3; Omieback = O; NOmieback = NO; NO2mieback = NO2;
        %% Second stage
        % Estimate reaction rates by multiplying rate coeffcients by Backward Euler concentrations
%         R1 = k1*NOmieback(1)*O3mieback(1);
%         R2 = J*NO2mieback(1);
%         R3 = k3*Omieback(1);
        %% Third stage
        % Time loop
        for i = 2:nsteps
            disp(['i = ' num2str(i)])
            % Iteration loop
            m = 0;
            np = 0;
            while np <= Np
                %disp(['np = ' num2str(np)])
                m = m + 1;
                P(i,1) = k1*NOmieback*O3mieback; % NO2
                P(i,2) = J*NO2mieback; % NO
                P(i,3) = J*NO2mieback; % O
                P(i,4) = k3*Omieback; % O3
                L(i,1) = J; % NO2, L - implicit loss coefficient (lambda)
                L(i,2) = k1*O3mieback; % NO
                L(i,3) = k3; % O
                L(i,4) = k1*NOmieback; % O3
                NO2mieback = (NO2mie(i-1) + dt*P(i,1))/(1 + dt*L(i,1));
                NO2mieforw = NO2mie(i-1) + dt*(P(i,1) - L(i,1));
                % Check convergence
                if any(NO2mieforw < 0) % || NOmieforw < 0 || Omieforw < 0 || O3mieforw < 0)
                    np = 0;
                else
                    np = np +1;
                end
                %% Six stage
                % Limit current backward Euler concentrations
                NO2mieback = min(NO2mieback, NO2miemax);
                % Update maximum backward Euler concentrations
                NO2miemax = max(NO2mieback, NO2mie(i-1));
            end
            %% Seventh stage
            NO2mie(i) = NO2mieback;
        end
end
%% Plot results
switch scheme
    case 1
        figure;
        end1 = 20;
        subplot(2,2,1); plot(O3forw(1:end1),'b'); title('O3')
        subplot(2,2,2); plot(Oforw(1:end1),'b'); title('O')
        subplot(2,2,3); plot(NOforw(1:end1),'b'); title('NO')
        subplot(2,2,4); plot(NO2forw(1:end1),'b'); title('NO2')
        %subplot(2,2,1); plot(t0:dt:nt,O3forw(1:end-1),'b'); title('O3')
        %subplot(2,2,2); plot(t0:dt:nt,Oforw(1:end-1),'b'); title('O')
        %subplot(2,2,3); plot(t0:dt:nt,NOforw(1:end-1),'b'); title('NO')
        %subplot(2,2,4); plot(t0:dt:nt,NO2forw(1:end-1),'b'); title('NO2')
    case 2
        figure;
        subplot(2,2,1); plot(O3back(1:end-1),'b'); title('O3')
        subplot(2,2,2); plot(Oback(1:end-1),'b'); title('O')
        subplot(2,2,3); plot(NOback(1:end-1),'b'); title('NO')
        subplot(2,2,4); plot(NO2back(1:end-1),'b'); title('NO2')
end