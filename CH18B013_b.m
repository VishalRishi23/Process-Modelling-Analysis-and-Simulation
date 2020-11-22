%fsolve
options=optimoptions('fsolve','functionTolerance',1e-30,'StepTolerance',1e-100,'OptimalityTolerance',1e-100);
Temp=fsolve(@equil,[500],options);

%Display the results
disp('The final temperature (Kelvin)')
disp(Temp)

%Function to evaluate temperature at equilibrium
function f = equil(T)
    global n0 n3 P;
    n0=100;
    n3=220;
    P=5000;
    z=0.8*n0;
    P_CO=((n0-z)/(n0+n3-(2*z)))*P;
    P_H2=((n3-(2*z))/(n0+n3-(2*z)))*P;
    P_M=((z)/(n0+n3-(2*z)))*P;
    Kp=P_M/(P_CO*(P_H2^2));
    f=Kp-((1.39E-4)*exp((21.225)+(9143.6/T)-(7.492*log(T))+((4.076E-3)*T)-((7.161E-8)*(T^2))));
end