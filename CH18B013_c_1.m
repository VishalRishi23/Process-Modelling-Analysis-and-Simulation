global n0 n3 n1 n2 Kp T P;
n1=100
n2=200
P=2000;
T=500;
Kp=(1.39E-4)*exp((21.225)+(9143.6/T)-(7.492*log(T))+((4.076E-3)*T)-((7.161E-8)*(T^2)));
n4=0
n5=0
n4vec=[];
n5vec=[];
%Successive substitution
for i = 1:100
    n0=n1+n4;
    n3=n2+n5;
    n4vec(i)=n4;
    n5vec(i)=n5;
    %fsolve
    options=optimoptions('fsolve','functionTolerance',1e-15,'StepTolerance',1e-15,'OptimalityTolerance',1e-15);
    zeta=fsolve(@equil,[50],options);
    n4=n0-zeta;
    n5=n3-(2*zeta);
    n6=zeta;
end

%Plot the figures
figure('Name','n4 vs iterations')
plot(n4vec,'-s')
xlabel('iterations')
ylabel('n4')
figure('Name','n5 vs iterations')
plot(n5vec,'-o')
xlabel('iterations')
ylabel('n5')

%Display the molar flow rates
disp('Value of n0 (kmol/h)')
disp(n0)
disp('Value of n1 (kmol/h)')
disp(n1)
disp('Value of n2 (kmol/h)')
disp(n2)
disp('Value of n3 (kmol/h)')
disp(n3)
disp('Value of n4 (kmol/h)')
disp(n4)
disp('Value of n5 (kmol/h)')
disp(n5)
disp('Value of n6 (kmol/h)')
disp(n6)
disp('Value of extent of reaction')
disp(zeta)

%Function to evaluate zeta at equilibrium
function f = equil(z)
    global n0 n3 Kp P T;
    T=500; 
    Kp=(1.39E-4)*exp((21.225)+(9143.6/T)-(7.492*log(T))+((4.076E-3)*T)-((7.161E-8)*(T^2)));
    P=2000;
    P_CO=((n0-z)/(n0+n3-(2*z)))*P;
    P_H2=((n3-(2*z))/(n0+n3-(2*z)))*P;
    P_M=((z)/(n0+n3-(2*z)))*P;
    f=(P_M/(P_CO*(P_H2^2)))-Kp;
end
    