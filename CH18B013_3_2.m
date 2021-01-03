tfinal = 200;
%Parameters
rho = 1000; rhoj = 800; cp = 2500; cpj = 5000; UA = 20000;
V = 0.25; Vj = 0.04;
Q = 25e-3; Qj = 5e-3; Tin = 350; Tjin = 298; cAin = 4000;
k0 = 800; EbyR = 4500; deltaH = -250000;
tau = V / Q; beta = UA / (Qj * rhoj * cpj);
pars = [rho rhoj cp cpj UA V Vj Q Qj Tin Tjin cAin k0 EbyR deltaH];
T = 270:0.1:748;
k = k0*exp(-EbyR./T);

% Graphical Method 
f2 = -deltaH*(k*tau./(1+k*tau))*Q*cAin + Q*rho*cp*(Tin-T)-(UA/(1+beta))*(T-Tjin);
f3 = f2./(Q*rho*cp*(Tin-Tjin)); % normalized
figure('units','normalized','outerposition',[0 0 1 1])
plot(T,f2,'LineWidth',2); refline(0,0); grid on;
xlabel('Tank Temperature (K)','FontSize',14); ylabel('f_2(T)','FontSize',14);

% Initial Guess is chosen from the above graphical method
[a, b] = find((abs(f3)<0.0007));
Tss0 = T(b);
Tjss0 = ((Qj*rhoj*cpj)*Tjin + UA*Tss0)/(Qj*cpj*rhoj + UA);
cAss0 = (Q/V)*cAin./(Q/V + k0*exp(-EbyR./Tss0));

% Use fsolve & plot in f2 vs. T plot
for i = 1:3
    sssolution = fsolve(@(y) cstrss(y,pars), [cAss0(i) Tss0(i) Tjss0(i)]);
    cAss(i) = sssolution(1); Tss(i) = sssolution(2); Tjss(i) = sssolution(3);
end
hold on
plot(Tss(1),0,'bo','MarkerSize',8,'MarkerFaceColor','blue')
hold on
plot(Tss(2),0,'bo','MarkerSize',8,'MarkerFaceColor','red')
hold on
plot(Tss(3),0,'bo','MarkerSize',8,'MarkerFaceColor','blue')

A=cAss;
B=Tss;
C=Tjss;

% Phase Plane plot for the Non-Linear Equations
figure('Name','Phase Plane Plot','units','normalized','outerposition',[0 0 1 1])
plot(cAss(1),Tss(1),'o','MarkerSize',18,'MarkerFaceColor','blue')
hold on
plot(cAss(2),Tss(2),'o','MarkerSize',18,'MarkerFaceColor','red')
hold on
plot(cAss(3),Tss(3),'o','MarkerSize',18,'MarkerFaceColor','blue')
hold on
cA0vec = [500:500:4000]; T0vec = [300:25:500]; Tj0 = 300;
for i = 1:numel(cA0vec)
    for j = 1:numel(T0vec)
        cA0 = cA0vec(i);
        T0 = T0vec(j);
        [t, cT] = ode15s(@(t,cT) cstrss(cT,pars), [0 tfinal], [cA0 T0 Tj0]);
        cAvec = cT(:,1); Tvec = cT(:,2);
        plot(cAvec, Tvec, 'o-')
        xlabel('Concentration of A (mol/m3)','FontSize',14); ylabel('Tank Temperature (K)','FontSize',14);
        hold on;
        plot(cA0,T0,'o','MarkerSize',8,'MarkerEdgeColor','green','MarkerFaceColor','green');
        hold on;
    end
end

% f2 vs. T plots for Tjin = 484K
Tjinvec = [484];
figure('Name','f2 vs. T for different Tjin','units','normalized','outerposition',[0 0 1 1])
for i = 1:numel(Tjinvec)
    Tjin = Tjinvec(i);
    f2 = -deltaH*(k*tau./(1+k*tau))*Q*cAin + Q*rho*cp*(Tin-T)-(UA/(1+beta))*(T-Tjin);
    f3 = f2/max(f2);
    plot(T,f2,'LineWidth',2); refline(0,0); grid on;
    hold on
    xlabel('Tank Temperature (K)','FontSize',14); ylabel('f_2(T)','FontSize',14); title(num2str(Tjin),'FontSize',14);
    [a, b] = find((abs(f3)<0.01));
    Tss0 = T(b);
    Tjss0 = ((Qj*rhoj*cpj)*Tjin + UA*Tss0)/(Qj*cpj*rhoj + UA);
    cAss0 = (Q/V)*cAin./(Q/V + k0*exp(-EbyR./Tss0));
    for j = 1:numel(b)
        pars = [rho rhoj cp cpj UA V Vj Q Qj Tin Tjin cAin k0 EbyR deltaH];
        sssolution = fsolve(@(y) cstrss(y,pars), [cAss0(j) Tss0(j) Tjss0(j)]);
        cAss(j) = sssolution(1); Tss(j) = sssolution(2); Tjss(j) = sssolution(3);
    end
end

% Bifurcation plot with Tj,in as the Bifurcation Parameter
Tin = 350;
clear Tjss 
T = 273:0.1:723;
k = k0*exp(-EbyR./T);
Tjinvec = [50 :3 :500];
Tjssmat = nan*ones(numel(Tjinvec),3);
for i = 1:numel(Tjinvec)
    Tjin = Tjinvec(i);
    f2 = -deltaH*(k*tau./(1+k*tau))*Q*cAin + Q*rho*cp*(Tin-T)-(UA/(1+beta))*(T-Tjin);
    f3 = f2/max(f2);
    b = find(abs(diff(sign(f2)))>0);
    Tss0 = T(b);
    Tjss0 = ((Qj*rhoj*cpj)*Tjin + UA*Tss0)/(Qj*cpj*rhoj + UA);
    cAss0 = (Q/V)*cAin./(Q/V + k0*exp(-EbyR./Tss0));
    for j = 1:numel(b)
        pars = [rho rhoj cp cpj UA V Vj Q Qj Tin Tjin cAin k0 EbyR deltaH];
        sssolution = fsolve(@(y) cstrss(y,pars), [cAss0(j) Tss0(j) Tjss0(j)]);
        cAss(i,j) = sssolution(1); Tss(i,j) = sssolution(2); Tjssmat(i,j) = sssolution(3);
    end
    Tjssmat;
    clear b pos
end   
figure('Name','Bifurcation Plot','units','normalized','outerposition',[0 0 1 1])
plot(Tjinvec, Tjssmat,'ro','MarkerSize',8,'MarkerFaceColor','red')
xlabel('Inlet Jacket Temperature (K)','FontSize',14); ylabel('Steady State Jacket Temperature (K)','FontSize',14); title('Bifurcation Plot','FontSize',14);

disp('Steady state 1 (cA,T,Tj):')
disp(A(1))
disp(B(1))
disp(C(1))
disp('Steady state 2 (cA,T,Tj):')
disp(A(2))
disp(B(2))
disp(C(2))
disp('Steady state 3 (cA,T,Tj):')
disp(A(3))
disp(B(3))
disp(C(3))

%Function definition
function f = cstrss(y,pars)
    cAss = y(1); Tss = y(2); Tjss = y(3);
    parsvec = num2cell(pars);
    [rho, rhoj, cp, cpj, UA, V, Vj, Q, Qj, Tin, Tjin, cAin, k0, EbyR, deltaH] = deal(parsvec{:});
    rA = k0*exp(-EbyR./Tss)*cAss;
    f(1,1) = ((Q/V)*(cAin-cAss)-rA);
    f(2,1) = (Q*rho*cp*(Tin - Tss) - deltaH*rA*V - UA*(Tss-Tjss))/(V*rho*cp);
    f(3,1) = (Qj*rhoj*cpj*(Tjin-Tjss) + UA*(Tss-Tjss))/(Vj*rhoj*cpj);
end