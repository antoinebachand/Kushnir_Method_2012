
%Normalized gas parameters
para.Ti = 1.15;
para.gamma = 1.4;
para.Z0 = 1;
para.ZT0 = 0;
para.R = (para.gamma-1)*para.Z0;
para.U = 0; %(para.gamma-1)*ZT0*T0;

% Heat exchange parameters
para.TRw = 1;
para.qr=10;
para.F0 = 4e-4;
para.Bi = 200;

% Leak parameters
para.kappa = 0.1;
para.FR = 15;
para.Pe = 1; 

% Injection scenario for one cycle
% In this exemple, 
%       injection is mr (Fie(1) = 1) from t = 0 to 2/12 (td(1) = 2)
%       injection is 0.2*mr (Fie(2) = 0.2) from t = 2/12 to 5/12 (td(2) = 5/12)
%       injection is 0 (Fie(3) = 0) from t = 5/12 to t = 1 (td(3) = 1)
%
% Important: length(para.Fie) = length(para.td)

para.mr = 0.9;  
para.Fie = [1 0.2 0];
para.td = [2 5 12]/12;
                         
SolverOptions.N = 20;
SolverOptions.beta = 1.05;
SolverOptions.Rp = 0;
SolverOptions.eps = 0.005;

NbCycles = 10;
[t,rho, Tair, P, Trock,rrock,ml,Nc] = fcnSolveKushnirNum(NbCycles,para,SolverOptions);

% perfect conductor kR -> \infty => F0 -> \infty and Bi -> 0
para.F0 = 1e6; 
para.Bi = 0;
[tinf,rhoinf, Tairinf, Pinf, Trockinf,rrockinf,mlinf,Ncinf] = fcnSolveKushnirNum(NbCycles,para,SolverOptions);

dt = 1/(12*60); 
% Zhou solution
% N.B. The function "fcnZhouSolution" is not well coded! It's only
% convinient for this particular injection scenario
[tZhou, rhoZhou, TZhou, PZhou, mlZhou,cnt] = fcnZhouSolution(NbCycles,para,dt);

tscale = 0.5;
figure(1)
subplot(2,1,1);
strcnt = int2str(cnt);
plot(tZhou*tscale,PZhou,':b',t*tscale,P,'--r','LineWidth',2);
hold on
plot(tinf*tscale,Pinf,'-k','LineWidth',0.2);
hold off
legend(strcat('Zhou (',strcnt,' iterations)'),...
            'Modified Kushnir (k_R = 4 W/mK)', ...
            'Modified Kushnir (k_R = \infty)','Location','southeast');
xlabel('t (days)');
ylabel('P/P_0');
%ylim([0 7]);
subplot(2,1,2)
plot(tZhou*tscale,mlZhou,':b',t*tscale,ml,'--r','LineWidth',2);
hold on
plot(tinf*tscale,mlinf,'k','LineWidth',0.2);
hold off
legend(strcat('Zhou (',strcnt,' iterations)'),...
        'Modified Kushnir (k_R = 4 W/mK)',...
        'Modified Kushnir (k_R = \infty)','Location','southeast');
xlabel('t (days)');
ylabel('Normalized leakage rate (m_l/m_c)');

print ('-f1', '-dpng', strcat('MultiCycles10daysPressureAndLeakageRate.png'));


figure(2)
subplot(2,1,1);
plot(tZhou*tscale,TZhou,':b',t*tscale,Tair,'--r','LineWidth',2);
hold on
plot(tinf*tscale,Tairinf,'k','LineWidth',0.5);
hold off
legend(strcat('Zhou (',strcnt,' iterations)'),...
        'Modified Kushnir (k_R = 4 W/mK)',...
          'Modified Kushnir (k_R = \infty)','Location','southeast');
xlabel('t (days)');
ylabel('T/T_0');
subplot(2,1,2);
plot(tZhou*tscale,rhoZhou,':b',t*tscale,rho,'--r',tinf*tscale,rhoinf,'k');
hold on
plot(tinf*tscale,rhoinf,'k');
hold off
legend(strcat('Zhou (',strcnt,' iterations)'),...
            'Modified Kushnir (k_R = 4 W/mK)',...
            'Modified Kushnir (k_R = \infty)','Location','southeast');
xlabel('t (days)');
ylabel('\rho/\rho_0');

print ('-f2', '-dpng', strcat('MultiCycles10daysTemperaturAndDensity.png'));




