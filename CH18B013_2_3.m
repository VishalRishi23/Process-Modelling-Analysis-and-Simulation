global Af D k;
Af=10; %Concentration of A in feed (mol/lit)
D=4/7; %Dilution rate (1/s)
k=[(5/6),(5/3),(1/6)]; %Rate constants
Cs=fsolve(@concentration,[0;0;0;0]);

disp('Steady state concentration of A (mol/lit)')
disp(Cs(1))
disp('Steady state concentration of B (mol/lit)')
disp(Cs(2))
disp('Steady state concentration of C (mol/lit)')
disp(Cs(3))
disp('Steady state concentration of D (mol/lit)')
disp(Cs(4))

J=jacobian(Cs);
[e_vectors,e_values]=eig(J);

disp('Eigen vlaues')
disp(e_values)
disp('Eigen vectors')
disp(e_vectors)

%Phase plane plot between Ca and Cb
figure
plot(Cs(1),Cs(2),'o','MarkerSize',10)
title('Phase plane plot between Cb and Ca')
xlabel('Ca (mol/lit)')
ylabel('Cb (mol/lit)')
for i=1:10
    hold on
    Ca=1+((4)*rand(1,1));
    Cb=0.3+((2.7)*rand(1,1));
    plot(Ca,Cb,'*','MarkerSize',20)
    [t,c]=ode45(@diff_eq,[0;200],[Ca;Cb]);
    plot(c(:,1),c(:,2),'-')
end
hold off

%Function definitions
function f=concentration(C)
    global Af D k;
    %Conservation equations
    f1=(D*(Af-C(1)))-(k(1)*C(1))-(2*k(3)*(C(1)^2));
    f2=(k(1)*C(1))-(k(2)*C(2))-(D*C(2));
    f3=(k(2)*C(2))-(D*C(3));
    f4=(k(3)*(C(1)^2))-(D*C(4));
    f=[f1;f2;f3;f4];
end

function J=jacobian(C)
    global D k;
    a=-(D+k(1)+(4*k(3)*C(1)));
    b=-(D+k(2));
    J=[a,0;k(1),b];
end

function dCdt=diff_eq(t,C)
    global Af D k;
    dC1=(D*(Af-C(1)))-(k(1)*C(1))-(2*k(3)*(C(1)^2));
    dC2=(k(1)*C(1))-(k(2)*C(2))-(D*C(2));
    dCdt=[dC1;dC2];
end

%Conservation equations:
%dA/dt = (Dil*(Af-A))-(k1*A)-(2*k3*(A^2))
%dB/dt = (k1*A)-(k2*B)-(Dil*B)
%dC/dt = (k2*B)-(Dil*C)
%dD/dt = (k3*(A^2))-(Dil*D)

%Stability:
%The eigen values of the Jacobian at the steady state
%are negative which indicates a stable steady state.

%Explanation:
%From the linear stability analysis, the eigen values of the Jacobian at the steady state 
%are negative which indicates a stable steady state. Hence, in the phase
%plane plot, whatever be the initial concentrations of A and B (indicated by the '*'), they
%finally reach their steady state concentrations, indicated by the 'o'.
