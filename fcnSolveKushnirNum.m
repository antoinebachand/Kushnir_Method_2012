function [t,rho, Tair, P, Trock,r,ml,Nc] = fcnSolveKushnirNum(NbCycles,para,SolverOptions)
%
%   This function will solve the differential coupled system of equations in Kushnir article
%   (2012 Temperature and pressure variations within compressed air energy
%   sorage cavern). The original equations are modified to include a
%   leakage rate. 
%
%   In this code, the parameter "mr" is defined as mr = mc*tp/rho_0/V 
%   where "tp" is the period of the CAES plant cycle. In Kushnir's article,
%   mr is defined as mr =  mc*t1/rho_0/V, where t1 represents the charging time.  
%
%   The code solve the differential system of equations for "NbCycles" identical  
%   CAES plant cycle (NbCycles HAS TO BE AN INTEGER). A CAES plant cycle
%   starts at t=0 and finiches at t = 1 (t is a normalized time)
%
%   The structure "para" is used to pass all the normalized parameters 
%   (eq. 17 of Kushnir) 
%
%       %Normalized gas parameters
%       para.Ti     -> Ti/T0
%       para.gama   -> cp0/cv0
%       para.R      -> (gamma-1)*Z0
%       para.U      -> (gamma-1)*ZTO*TO
%       
%       %Normalized heat exchange parameters
%       para.TRw    -> TRw0/T0 %(Normalized cavern wall temperature)
%       para.qr     -> hc*Ac/cv0/mc
%       para.F0     -> alpha_R*tp/Rw/Rw
%       para.Bi     -> hc*Rw/kr
%
%       % Normalized leak parameters (eq 8 of Zhou); 
%       para.kappa   -> 1.0967e-2*kr*H*(para.gamma-1)*cv0*para.rho0*P0/mc/mu;
%       para.FR      -> kr*tp*P0/mu/phi/Rw^2;
%       para.Pe      -> Pe/P0
%
%       % Injection scenario for one cycle
%       
%       para.Fie     -> Fi or Fe; Fie > 0 (injection) Fie<0 (Discharge) Fie = 0 (Storage)
%                        Ex : Fie = x => injection rate is x*mc
%       para.mr      -> mc*tp/rho0/V where mc is the mass injection rate 
%                        (mr is different from Kushnir!!!)
%       td indicates the normalized time where Fie is changing
%           Exemple: 
%                 If td = [2 5 12]/12 and Fie = [1 0.2 0]
%                 then:
%                   From t=0 to t=2/12 => Fie = 1
%                   From t=2/12 to t=5/12 => Fie = 0.2
%                   From t=5/12 to t=1 (end of cycle) => Fie = 0
% 
%   The structure "SolverOptions" is use to pass the numerical grid info 
%   for the transfromed uniform grid in the rock (eq. 18 of Kushnir)
%       
%   SolverOptions.N     -> N is the number of cells for the 
%                           discretized Temperature in the rock
%   SolverOptions.beta  -> beta is the parameter controlling the uniform
%                           grid transformation (eq. 18) (1 < beta < \infty)
%   SolverOptions.Rp    ->  Rp is the normalized penetration radius
%                           if Rp = 0 then Rp is estimated with eq A.7 of
%                           Kushnir
%   SolverOptions.eps   -> epsilon in eq. A7
%
%   The coupled system of equations is store in variable Y:
%   Y[1] -> Temperature in the rock at r = 1 (j=1) (at cavern wall)
%   Y[2] -> Temperature in the rock at r = r (j=2) (eq 24 of Kushnir; j starts at 1 in this code!)
%     :
%   Y[N+1] -> Temperature in the rock at r = Rp
%   Y[N+2] -> air temperature (eq. 27)
%   Y[N+3] -> density (Normalized version of eq 1 wher we substract a  
%                       normalized leakage rate = ml/mc (with ml eq 8 of
%                       Zhou)
%
%   THE FUNCTION RETURNS:
%
%   t   : the time where the solution was evaluated   
%   rho : the normalized air density at different time
%   Tair: the normalized air temperature
%   Trock: the rock temperature at r (N+1 values of r: r=1..r=Rp)
%   ml: normalized leakage rate
%   Nc: index to indicate the end of cycle 
%           Nc(1) -> t(Nc(1)) = 1
%           Nc(2) -> t(Nc(2)) = 2
%             :
%           Nc(NbCycles) -> t(Nc(NbCycles)) = NbCycles 
%
%
% ======================================================================
% Copyright (c) November 2021, Bernard Doyon (bdoyon@cegepgarneau.ca)
% ======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Normalized gas parameters
    Ti = para.Ti;
    gamma =para.gamma;
    R = para.R;
    U = para.U;
    
    % Heat exchange parameters
    TRw = para.TRw; 
    qr = para.qr;
    F0 = para.F0;
    Bi = para.Bi;
    
    % Leak parameters
    kappa = para.kappa;
    FR = para.FR;
    Pe = para.Pe;
    
    % Injection scenario for one cycle
    Fie = para.Fie;
    td = para.td;
    mr = para.mr;   
            
    % Numerical options for the grid size  
    N = SolverOptions.N;    
    beta = SolverOptions.beta;
    flagRp = floor(SolverOptions.Rp);
    Rp = SolverOptions.Rp;
    eps = SolverOptions.eps;

    % initial condition
    Yinitial = ones(1,N+1)*TRw;  % Y(1:N+1): Rock temperature at N+1 values of r 
    Yinitial = [Yinitial 1 1];   % Y(N+2):   Air temp
                                 % Y(N+3):   rho
     
    eta = linspace(0,1,N+1);    % eta parameter from 0 to 1 (N+1 values)
    eta_p = zeros(1,N+1);       % first derivative of eta with respect to r
    eta_pp = zeros(1,N+1);      % second derivative of eta with respect to r
    r = zeros(NbCycles,N+1);    % the r values in the rock whete 
                                %the temperatur is evaluated. From cycle 
                                % to cycle the "r" values are different 
                                %since Rp is different                              
    Nd = length(para.td);
    Nc = zeros(1,NbCycles);
    
    % main loop
    for n = 1:NbCycles
        % Evaluate Rp from eq A.7 of Kushnir
        if not(flagRp)
            xi_n = sqrt(2*pi/(F0*n));
            Rp = 1+sqrt(2)/xi_n*log(1/eps);
        end
        r(n,:) = fcnrj(eta,beta,Rp);
        for j = 1:N+1
            eta_p(j) = fcneta_p(r(n,j),beta,Rp);  % calculate first derivative of eta
            eta_pp(j) = fcneta_pp(r(n,j),beta,Rp); %  calculate second derivative of eta
        end     
        for i=1:Nd % for every discontinuous points of Fi or Fe
            if i ==1 
                tspan =[0 td(i)]; % the time interval where the solution is calculated
            else
                tspan = [td(i-1) td(i)];
            end
            if (Fie(i) >0) %charge -> use the value of gamma
                % stiff solver 
                [tspan,Yspan] = ode15s(@(t,Y) odefcnChargeDischargeStorage(t,Y,n,Fie(i),mr,Ti,gamma,R,U,F0,qr,Bi,kappa,FR,Pe,N,r(n,:),eta_p,eta_pp),tspan,Yinitial);
            else % Discharge or Storage - > set gamma to zero
                [tspan,Yspan] = ode15s(@(t,Y) odefcnChargeDischargeStorage(t,Y,n,Fie(i),mr,Ti,0,R,U,F0,qr,Bi,kappa,FR,Pe,N,r(n,:),eta_p,eta_pp),tspan,Yinitial);
            end
            Yinitial = Yspan(end,:); % set the inttial condition for next time interval 
            % store the calculated values
            % the number of calculated values will change from cycle to
            % cycle
            if ((i == 1)&&(n==1))
                t = (n-1)+tspan;
                rho = Yspan(:,N+3);
                P = Yspan(:,N+3).*Yspan(:,N+2);
                Tair = Yspan(:,N+2);
                Trock = Yspan(:,1:N+1);
            else
                t = [t;(n-1)+tspan]; %#ok<*AGROW>
                rho = [rho;Yspan(:,N+3);];
                P = [P;Yspan(:,N+3).*Yspan(:,N+2)];
                Tair = [Tair;Yspan(:,N+2)];
                Trock = [Trock;Yspan(:,1:N+1)];
            end
            Nc(n) = Nc(n) + length(tspan);
        end     
    end
    % evaluate the leak for all t values in the cycle 
    ml = fcnCalculateLeak(t,kappa,FR,Tair,P,Pe);
end



function dYdt = odefcnChargeDischargeStorage(t,Y,n,Fie,mr,Ti,gamma,R,U,F0,qr,Bi,kappa,FR,Pe,N,r,eta_p,eta_pp)
   
    % Charge: Fie >0
    % Discharge Fie <0 AND gamma = 0
    % Storage Fie = 0
    dYdt = zeros(N+3,1);
    
    P = Y(N+3)*Y(N+2);
    ml = fcnCalculateLeak(n-1+t,kappa,FR,Y(N+2),P,Pe);
    
    dYdt(1) = F0*(2*N^2*eta_p(1)^2*(Y(2)-Y(1))+ Bi*(eta_pp(1)/eta_p(1) + 1/r(1) - 2*N*eta_p(1))*(Y(1)-Y(N+2)));
    for j=2:N
       dYdt(j) = F0*(N^2*eta_p(j)^2*(Y(j+1)-2*Y(j)+Y(j-1)) + N*(eta_pp(j)+eta_p(j)/r(j))*(Y(j+1)-Y(j-1))/2); 
    end
    dYdt(N+1) = F0*(2*N^2*eta_p(N+1)^2*(Y(N)-Y(N+1)));
    dYdt(N+2) = (Fie-ml)*mr/Y(N+3)*(gamma*Ti+(R-gamma)*Y(N+2)+U*Y(N+3)) + mr/Y(N+3)*qr*(Y(1)-Y(N+2));
    dYdt(N+3) = (Fie-ml)*mr;
end

function [ml] = fcnCalculateLeak(t,kappa,FR,T,P,Pe)
    % 
    % NOT SURE ABOUT THIS VALUE!
    %
    % figure 5 of Zhou shows a constant air leakage for the intial
    % iteration (after 6 hours). If Re depends on time, the air leakage 
    % should not be constant!
    ReRw = 1 + sqrt(FR*P.*t);
%    ReRw = 1 + sqrt(FR*(P+1)/2);
%    Pmean = mean(P);
 %   ReRw = 1 + sqrt(FR*Pmean);
    ml = kappa*(P.^2-Pe.^2)./T./log(ReRw);
    
    ml(isnan(ml))=0;
    ml(isinf(ml))=0;
end


function rj = fcnrj(eta,beta,Rp)
    % the values of r for the uniform interval of eta (0<eta<1) 
    rj = 1+(Rp-1)*(beta+1-(beta-1)*((beta+1)/(beta-1)).^(1-eta))./(1+((beta+1)/(beta-1)).^(1-eta));
end

function eta_p = fcneta_p(r,beta,Rp)
    % first derivative of eta with respect to r
  eta_p = -(-1/(Rp-1)/(beta-1+(r-1)/(Rp-1))-(beta+1-(r-1)/(Rp-1))/(beta-1+(r-1)/(Rp-1))^2/(Rp-1))/(beta+1-(r-1)/(Rp-1))*(beta-1+(r-1)/(Rp-1))/log(((beta+1)/(beta-1)));
end

function eta_pp = fcneta_pp(r,beta,Rp)
    % second derivative of eta with respect to r
   eta_pp = -(2 / (Rp - 1) ^ 2 / (beta - 1 + (r - 1) / (Rp - 1)) ^ 2 + 2 * (beta + 1 - (r - 1) / (Rp - 1)) / (beta - 1 + (r - 1) / (Rp - 1)) ^ 3 / (Rp - 1) ^ 2) / (beta + 1 - (r - 1) / (Rp - 1)) * (beta - 1 + (r - 1) / (Rp - 1)) / log(((beta + 1) / (beta - 1))) - (-1 / (Rp - 1) / (beta - 1 + (r - 1) / (Rp - 1)) - (beta + 1 - (r - 1) / (Rp - 1)) / (beta - 1 + (r - 1) / (Rp - 1)) ^ 2 / (Rp - 1)) / ((beta + 1 - (r - 1) / (Rp - 1)) ^ 2) * (beta - 1 + (r - 1) / (Rp - 1)) / log(((beta + 1) / (beta - 1))) / (Rp - 1) - (-1 / (Rp - 1) / (beta - 1 + (r - 1) / (Rp - 1)) - (beta + 1 - (r - 1) / (Rp - 1)) / (beta - 1 + (r - 1) / (Rp - 1)) ^ 2 / (Rp - 1)) / (beta + 1 - (r - 1) / (Rp - 1)) / (Rp - 1) / log(((beta + 1) / (beta - 1)));
end