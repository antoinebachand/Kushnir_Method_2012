%Caverne parameters
Rw = [2.8068 20 2.5]; 
Ac = [211.5 25000 15.708]; 
V = [222.75 141000 19.635];
H = [9 112.2 1];

n = 1; % n = 1 => Kamioka CAES
       % n = 2 => Huntorf plant
       % n = 3 => shut in test of Shou

Rw = Rw(n);
Ac = Ac(n);
V = V(n);
H = H(n);


% The air  parameters are for the Kamioka CAES 

% gas parameters

Ti = 30 + 273.15;   % Injected air temperature (K) 
T0 = 20 + 273.15;   % Initial air temperature (K)
P0 = 1.0133e5;      % Initial gaz pressure (Pa) 
cp0 = 1000;         % Constant pressure specific heat (J/Kg/K)
cv0 = 717;          % Constant volume specific heat (J/kg/K)

% Rock parameters

TRw0 = 20 + 273.15; % Initial rock temparature (K)
Pe = 1.0133e5;      % Air pressur at reservoir edge (Pa)
rhoR = 2100;        % Rock density (kg/m^3)
cpR = 840;          % Constant pressure specific heat of rock (J/kg/K)
hc = 30;            % Heat transfer coefficient  (W/m2/K)
kr = 5e-14;         % Calibrated permeability of rock (m^2)
kR = 4;             % thermal conductivity (W/m/K)
phi = 0.1;          % Porosity
mu = 1.79e-5;       % Air viscosity (Pa s)


% Injection scenario
% (mc is the maximal injection rate for one cycle)
mc = 18.38/60;          % Injection mass rate (kg/s) 
tp = 12*3600;           % Period of the cycle  (s)
td = [2.2 5.3 12]*3600; % Time where Fi is changing (s)
                        % (in this exemple Fi is changin at 2.2 hours and
                        % 5.3 hours and 12 hours (end of cycle)
     
Fie = [1 13.90/18.38 0]; % Normalized injection rate
                            % In this exemple, 
                            % From 0 to 2.2h -> Fi = 1 (meaning mc is 18.38 kg/s)
                            % From 2.2h to 5.3h -> Fi = 13.9/18.38 (meaning mc = 13.9 kg/s)
                            % From 5.2h to 12h  -> Fe = 0

% Initialize the structure "para" to pass to fcnSolveKushnirNum
para.rho0 = P0/(cp0-cv0)/T0;
para.Fie = Fie;
para.mr = mc*tp/para.rho0/V;        %mr is different from Kushnir
para.Ti = Ti/T0;
para.gamma = cp0/cv0;
para.qr=hc*Ac/cv0/mc;
para.kappa = 1.0967e-2*kr*H*(para.gamma-1)*cv0*para.rho0*P0/mc/mu;
para.FR = kr*tp*P0/mu/phi/Rw^2;
para.alphaR = kR/rhoR/cpR;
para.F0 = para.alphaR*tp/Rw^2;
para.Bi = hc*Rw/kR;
para.Z0 = 1;    % should be specifically adapted for hydrogen  
para.ZT0 = 0;   % should be specifically adapted for hydrogen 
para.R = (para.gamma-1)*para.Z0;
para.U = (para.gamma-1)*para.ZT0*T0;
para.td = td/tp;
para.T0 = T0;
para.TRw = TRw0/T0;
para.Pe = Pe/P0;

SolverOptions.N = 20;
SolverOptions.beta = 1.05;
SolverOptions.Rp = 0;
SolverOptions.eps = 0.005;

% Solve for one cycle with modified Kushnir equations
NbCycles = 1;
[t,rho, Tair, P, Trock,rrock,ml,Nc] = fcnSolveKushnirNum(NbCycles,para,SolverOptions);

% Compare with perfect conductor
kR = 1e6;
para.alphaR = kR/rhoR/cpR;
para.F0 = para.alphaR*tp/Rw^2;
para.Bi = hc*Rw/kR;
[tinf,rhoinf, Tairinf, Pinf, Trockinf,rrockinf,mlinf,Ncinf] = fcnSolveKushnirNum(NbCycles,para,SolverOptions);

%Zhou solution
dt = 1/(tp/60);
[tZhou, rhoZhou, TZhou, PZhou, mlZhou,cnt] = fcnZhouSolution(NbCycles,para,dt);


tscale = 12; % time scale for the x axis 
figure(1)
subplot(2,1,1);
plot(tZhou*tscale,PZhou,':b',t*tscale,P,'--r','LineWidth',1.5);
hold on
plot(tinf*tscale,Pinf,'k','LineWidth',0.5);
hold off
strcnt = int2str(cnt);
legend(strcat('Zhou (',strcnt,' iterations)'),'Modified Kushnir (k_R = 4 W/mK)', 'Modified Kushnir (k_R = \infty)');
xlabel('t (hours)');
ylabel('P/P_0');
%ylim([0 7]);
subplot(2,1,2)
plot(tZhou*tscale,mlZhou*mc,':b',t*tscale,ml*mc,'--r','LineWidth',1.5);
hold on
plot(tinf*12,mlinf*mc,'k','LineWidth',0.5);
hold off
legend(strcat('Zhou (',strcnt,' iterations)'),'Modified Kushnir (k_R = 4 W/mK)', 'Modified Kushnir (k_R = \infty)');
xlabel('t (hours)');
ylabel('Leakage rate (kg/s)');

print ('-f1', '-dpng', strcat('Figures6and7Zhou.png'));

figure(2)
subplot(2,1,1);
plot(tZhou*tscale,TZhou,':b',t*tscale,Tair,'--r',tinf*tscale,Tairinf,'.k');
legend(strcat('Zhou (',strcnt,' iterations)'),'Modified Kushnir (k_R = 4 W/mK)', 'Modified Kushnir (k_R = \infty)');
xlabel('t (hours)');
ylabel('T/T_0');
subplot(2,1,2);
plot(tZhou*tscale,rhoZhou,':b',t*tscale,rho,'--r',tinf*tscale,rhoinf,'.k');
legend(strcat('Zhou (',strcnt,' iterations)'),'Modified Kushnir (k_R = 4 W/mK)', 'Modified Kushnir (k_R = \infty)');
xlabel('t (hours)');
ylabel('\rho/\rho_0');

figure(3) % plotting the rock temparature ; r is the normalized radius (compare to the cavern radius)
subplot(2,1,1);
plot(t*tscale,Trock(:,1),t*tscale,Trock(:,4),t*tscale,Trock(:,8),t*tscale,Trock(:,12),t*tscale,Trock(:,15))
xlabel('t (hours)');
ylabel('T_{rock}/T_0');
rstr1 = num2str(rrock(1));
rstr2 = num2str(rrock(4));
rstr3 = num2str(rrock(8));
rstr4 = num2str(rrock(12));
rstr5 = num2str(rrock(15));
legend(strcat('r =',rstr1),strcat('r =',rstr2),strcat('r =',rstr3),strcat('r =',rstr4),strcat('r =',rstr5));
subplot(2,1,2);
plot(rrock,Trock(1,:),rrock,Trock(end,:));
xlabel('r/R_w');
ylabel('T_{rock}/T_0');
legend('t = 0', 't = 12h','Location','southeast');


