global n0 n3 P T Kp;
n0=100;
n3=220;
P=1000;
T=600;
Kp=(1.39E-4)*exp((21.225)+(9143.6/T)-(7.492*log(T))+((4.076E-3)*T)-((7.161E-8)*(T^2)));

%fsolve
options=optimoptions('fsolve','functionTolerance',1e-30,'StepTolerance',1e-100,'OptimalityTolerance',1e-100);
zeta=fsolve(@equil,[50],options);
n4=n0-zeta;
n5=n3-(2*zeta);
n6=zeta;

%Display the results
disp('Value of n6 (kmol/h)')
disp(n6)

%Function to evaluate zeta at equilibrium
function f = equil(z)
    global n0 n3 Kp P T;
    n0=100;
    n3=220;
    T=600;
    Kp=(1.39E-4)*exp((21.225)+(9143.6/T)-(7.492*log(T))+((4.076E-3)*T)-((7.161E-8)*(T^2)));
    P=1000;
    P_CO=((n0-z)/(n0+n3-(2*z)))*P;
    P_H2=((n3-(2*z))/(n0+n3-(2*z)))*P;
    P_M=((z)/(n0+n3-(2*z)))*P;
    f=(P_M/(P_CO*(P_H2^2)))-Kp;
end