global n1 n2 P T Kp;
n1=100;
n2=200;
P=5000;
T=500;
Kp=(1.39E-4)*exp((21.225)+(9143.6/T)-(7.492*log(T))+((4.076E-3)*T)-((7.161E-8)*(T^2)));

%fsolve
options=optimoptions('fsolve','functionTolerance',1e-30,'StepTolerance',1e-100,'OptimalityTolerance',1e-100,'MaxFunctionEvaluations',1000);
flow=fsolve(@flowsheet,[0;0;0;0;0;50],options);

%Display the results
disp('Value of n0 (kmol/h)')
disp(flow(1))
disp('Value of n1 (kmol/h)')
disp(n1)
disp('Value of n2 (kmol/h)')
disp(n2)
disp('Value of n3 (kmol/h)')
disp(flow(2))
disp('Value of n4 (kmol/h)')
disp(flow(3))
disp('Value of n5 (kmol/h)')
disp(flow(4))
disp('Value of n6 (kmol/h)')
disp(flow(5))
disp('Value of extent of reaction')
disp(flow(6))

%Function to evaluate molar flow rates
function f = flowsheet(n)
    global n1 n2 Kp P T;
    T=500;
    n1=100;
    n2=200;
    Kp=(1.39E-4)*exp((21.225)+(9143.6/T)-(7.492*log(T))+((4.076E-3)*T)-((7.161E-8)*(T^2)));
    P=5000;
    P_CO=((n(1)-n(6))/(n(1)+n(2)-(2*n(6))))*P;
    P_H2=((n(2)-(2*n(6)))/(n(1)+n(2)-(2*n(6))))*P;
    P_M=((n(6))/(n(1)+n(2)-(2*n(6))))*P;
    f(1)=n1+n(3)-n(1);
    f(2)=n2+n(4)-n(2);
    f(3)=n(3)+n(6)-n(1);
    f(4)=n(4)+(2*n(6))-n(2);
    f(5)=n(5)-n(6);
    f(6)=Kp-(P_M/(P_CO*(P_H2^2)));
end