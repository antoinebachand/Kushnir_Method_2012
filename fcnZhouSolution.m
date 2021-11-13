function [t_all,rho_all,T_all,P_all,ml_all,cnt] = fcnZhouSolution(NbCycles,para,dt)

    % Only injection and first storage is coded.
    % The injection phase has two injection rates:
        % 0 ... t1 = mr
        % t1 ... t1_2 = mr_2
        % t1_2 ... t2 = storage
   

    mr = para.mr;
    Fie = para.Fie;
    Ti = para.Ti;
    gamma =para.gamma;
    qr = para.qr;
    FR = para.FR;
    kappa = para.kappa;
    td = para.td;
    TRw = para.TRw;
    Pe = para.Pe;

    
    N = 1/dt;
    t = linspace(0,1,N);
    
    flag = 1;
    ml_k  = 1e8;
    cnt=0;
    cnt_max = 50;
    % first part of Zhou solution
    
    rho_ini = 1;
    T_ini = 1;

    for n = 1:NbCycles
        delta = 1;
        cnt=0;
        ml = zeros(1,N);
        while((delta>0.01)&&(cnt<cnt_max))
            [rho,rho_av] = odefcnAirDensity(t,td,Fie,mr,ml,N,rho_ini);
            [T,flag] = fcnTemp(t,td,mr,Fie,Ti,gamma,qr,TRw,[0 0 0],rho_av,T_ini);%XIA
            P = rho.*T;
            P(P<=0)=Pe;
            [ml] = fcnCalculateLeak((n-1)+t,kappa,FR,T,P,Pe);
            delta = abs((ml(end)-ml_k)/ml_k);
            ml_k = ml(end);
            cnt=cnt+1;
        end
        if (cnt == cnt_max)
            error('Increase number of itetaration for Zhou solution')
        end
        
        %Adjustment iteration for Zhou solution
        
        delta = 1;
        cntB=0;
        while((delta>0.01)&&(cntB<cnt_max))
            [rho,rho_av] = odefcnAirDensity(t,td,Fie,mr,ml,N,rho_ini);
            ml_av = fcnCalculateMeanLeak(ml,flag);
            [T,flag] = fcnTemp(t,td,mr,Fie,Ti,gamma,qr,TRw,ml_av,rho_av,T_ini);%XIA
            P = rho.*T;
            P(P<=0)=Pe;
            [ml] = fcnCalculateLeak((n-1)+t,kappa,FR,T,P,Pe);
            delta = abs((ml(end)-ml_k)/ml_k);
            ml_k = ml(end);
            cntB=cntB+1;
        end
        if (cntB == cnt_max)
            error('Increase number of itetaration for Zhou solution')
        end
        
        % initial condition for next cycle
        rho_ini = rho(end);
        T_ini = T(end);
        
        if n ==1
            rho_all = rho;
            T_all = T;
            t_all = t;
            ml_all = ml;
            P_all = P;
        else
            t_all = [t_all (n-1)+t]; %#ok<*AGROW>
            rho_all = [rho_all rho];
            T_all = [T_all T];
            ml_all = [ml_all ml];
            P_all = [P_all P];
        end
    end
    cnt = cnt+cntB;
end



function [rho,rho_av] = odefcnAirDensity(t,td,Fie,mr,ml,N,rho_ini)
    rho = zeros(1,N);
    sum_ml = 0;
    NbFi = length(td);
    
    sum_mr = 0;
    rho_av = 0;
    flag = 1;
    dt = t(2)-t(1);
    for n = 1:(NbFi-1) 
        for j = flag:N
            if t(j)<= td(n)
                sum_ml = sum_ml + ml(j)*dt;
                if (n == 1)
                    rho(j) = rho_ini + sum_mr + Fie(n)*mr*t(j) - mr*sum_ml;
                    
                else
                    rho(j) = rho_ini + sum_mr + Fie(n)*mr*(t(j)-td(n-1)) - mr*sum_ml;
                end
            else
                if n ==1
                    sum_mr = sum_mr + Fie(n)*mr*(td(n));
                    rho(rho<0) = 1;
                    rho_av(1) = mean(rho(1:j-1));
                else
                    sum_mr = sum_mr + Fie(n)*mr*(td(n)-td(n-1));
                    rho(rho<0) = 1;
                    rho_av(n) = mean(rho(flag:j-1));
                end
                flag = j;
                break;
            end
        end
    end
    for j = flag:N
        sum_ml = sum_ml + ml(j)*dt; 
        rho(j) = rho_ini + sum_mr - mr*sum_ml;
    end
    rho(rho<0) = 1;
    rho_av(n+1)= mean(rho(flag:N));
end


function [T,flag] = fcnTemp(t,td,mr,Fie,Ti,gamma,qr,TRw,ml_av,rho_av,T_ini)
    alpha1 = (Fie(1)*gamma*Ti + qr*TRw)/(ml_av(1)*(1-gamma)-Fie(1)-qr);
    alpha2 = (Fie(2)*gamma*Ti + qr*TRw)/(ml_av(2)*(1-gamma)-Fie(2)-qr);
    alpha3 = (qr*TRw)/(ml_av(3)*(1-gamma)-qr);
    beta1 = mr/rho_av(1)*(ml_av(1)*(1-gamma)-Fie(1)-qr);
    beta2 = mr/rho_av(2)*(ml_av(2)*(1-gamma)-Fie(2)-qr);
    beta3 = mr/rho_av(3)*(ml_av(3)*(1-gamma)-qr);

    T = zeros(1,length(t));
    for j = 1:length(t)
        if t(j)<= td(1)
            T(j) = (T_ini + alpha1)*exp(beta1*(t(j)))-alpha1;
            T1 = T(j);
            flag(1) = j;
        elseif t(j) <=td(2)
            T(j) = (T1 + alpha2)*exp(beta2*(t(j)-td(1)))-alpha2;
            T2 = T(j);
            flag(2) = j;
        else
            T(j) = (T2 + alpha3)*exp(beta3*(t(j)-td(2)))-alpha3;
        end
    end
end

function [ml] = fcnCalculateLeak(t,kappa,FR,T,P,Pe)
    ReRw = 1 + sqrt(FR*P.*t);
%    ReRw = 1 + sqrt(FR*(P+1)/2);
%    Pmean = mean(P);
 %   ReRw = 1 + sqrt(FR*Pmean);
    ml = kappa*(P.^2-Pe.^2)./T./log(ReRw);
    ml(isnan(ml)) = 0;
end

function [ml_av] = fcnCalculateMeanLeak(ml,flag)
    ml_av(1) = mean(ml(1:flag(1)));
    ml_av(2) = mean(ml(flag(1)+1:flag(2)));
    ml_av(3) = mean(ml(flag(2)+1:end));
end

