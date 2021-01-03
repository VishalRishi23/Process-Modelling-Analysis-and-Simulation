global Y b k1 a umax km x2f Dinit x3max;
%Parameters
Y=0.4;
b=5.5e-5;
k1=0.04545;
a=2.2;
umax=1.33e-4;
km=1.2;
x2f=20;
Dinit=5.611e-5;
x3max=50;
x_ss1_init=[0;x2f;0];
x_ss2_init=steadystate(Dinit);

%Phase plane plot between x1 and x2
figure('Name','Phase plane','units','normalized','outerposition',[0 0 1 1])
plot(x_ss2_init(1),x_ss2_init(2),'o','MarkerSize',18,'MarkerFaceColor','red')
title('Phase plane plot between x1 and x2')
xlabel('x1 (kg/m3)')
ylabel('x2 (kg/m3)')
for i=1:15
    hold on
    x1_init=0+((5)*rand(1,1));
    x2_init=5+((15)*rand(1,1));
    x3_init=0+((19)*rand(1,1));
    plot(x1_init,x2_init,'o','MarkerSize',18,'MarkerFaceColor','green')
    [t,c]=ode15s(@diff_eq,[0;1000000],[x1_init;x2_init;x3_init]);
    plot(c(:,1),c(:,2),'o-','LineWidth',2)
end
hold off

Drange=linspace(0.5,7.,100)/1e5;
e_val1=[];
e_val2=[];
for j=1:100
    e_val1=[e_val1;eigen_values(Drange(j),1)];
end
x_ss1=[zeros(1,100);20*ones(1,100);zeros(1,100)];
for j=1:100
    e_val2=[e_val2;eigen_values(Drange(j),2)];
end
x_ss2=cellfun(@steadystate,num2cell(Drange),'uni',0);
x_ss2=cell2mat(x_ss2);

%Plot of Eigen values vs the Dilution rate for the trivial steady state
figure('Name','Eigen values vs D','units','normalized','outerposition',[0 0 1 1])
subplot(3,1,1)
plot(Drange,e_val1(:,1),'LineWidth',2)
title('Plot of Eigen value 1 vs the Dilution rate for the trivial steady state')
xlabel('Dilution rate (/s)')
ylabel('Eigen value 1')
subplot(3,1,2)
plot(Drange,e_val1(:,2),'LineWidth',2)
title('Plot of Eigen value 2 vs the Dilution rate for the trivial steady state')
xlabel('Dilution rate (/s)')
ylabel('Eigen value 2')
subplot(3,1,3)
plot(Drange,e_val1(:,3),'LineWidth',2)
title('Plot of Eigen value 3 vs the Dilution rate for the trivial steady state')
xlabel('Dilution rate (/s)')
ylabel('Eigen value 3')

%Plot of Eigen values vs the Dilution rate for the non-trivial steady state
figure('Name','Eigen values vs D','units','normalized','outerposition',[0 0 1 1])
subplot(3,1,1)
plot(Drange,e_val2(:,1),'LineWidth',2)
title('Plot of Eigen value 1 vs the Dilution rate for the non-trivial steady state')
xlabel('Dilution rate (/s)')
ylabel('Eigen value 1')
subplot(3,1,2)
plot(Drange,e_val2(:,2),'LineWidth',2)
title('Plot of Eigen value 2 vs the Dilution rate for the non-trivial steady state')
xlabel('Dilution rate (/s)')
ylabel('Eigen value 2')
subplot(3,1,3)
plot(Drange,e_val2(:,3),'LineWidth',2)
title('Plot of Eigen value 3 vs the Dilution rate for the non-trivial steady state')
xlabel('Dilution rate (/s)')
ylabel('Eigen value 3')

%Bifurcation plot of the trivial steady state with Dilution rate as the
%bifurcation parameter
figure('Name','Bifurcation plot 1','units','normalized','outerposition',[0 0 1 1])
plot(Drange,x_ss1(1,:),'LineWidth',2)
hold on
plot(Drange,x_ss1(2,:),'LineWidth',2)
hold on 
plot(Drange,x_ss1(3,:),'LineWidth',2)
title('Bifurcation plot of the trivial steady state')
xlabel('Dilution rate (/s)')
ylabel('Steady State Concentration (kg/m3)')
legend('x1','x2','x3')
hold off

%Bifurcation plot of the non-trivial steady state with Dilution rate as the
%bifurcation parameter
figure('Name','Bifurcation plot 2','units','normalized','outerposition',[0 0 1 1])
plot(Drange,x_ss2(1,:),'LineWidth',2)
hold on
plot(Drange,x_ss2(2,:),'LineWidth',2)
hold on 
plot(Drange,x_ss2(3,:),'LineWidth',2)
title('Bifurcation plot of the non-trivial steady state')
xlabel('Dilution rate (/s)')
ylabel('Steady State Concentration (kg/m3)')
legend('x1','x2','x3')
hold off

disp('Trivial steady state for D = 5.611e-5 (x1,x2,x3):')
disp(x_ss1_init)
disp('Non-trivial steady state D = 5.611e-5 (x1,x2,x3):')
disp(x_ss2_init)

%Bifurcation analysis:
%For trivial steady state:
%Dilution rate <= 6.803e-5, Saddle state
%Dilution rate > 6.803e-5, Stable state
%For non-trivial steady state:
%Dilution rate <= 6.7374e-5, Saddle state
%Dilution rate > 6.7374e-5, Stable state

%Function definitions
function y=steadystate(D)
global Y b k1 a umax km x2f x3max;
options=optimoptions('fsolve','functionTolerance',1e-30,'StepTolerance',1e-100,'OptimalityTolerance',1e-100,'MaxFunctionEvaluations',1000);
y=fsolve(@solver,[0;0;0],options);
    function f=solver(x)
        f1=(Y*(x2f-x(2)))-x(1);
        f2=(((a*D)+b)*x(1))-(D*x(3));
        f3=((umax*(1-(x(3)/x3max))*x(2))/(km+x(2)+(k1*(x(2)^2))))-D;
        f=[f1;f2;f3];
    end
end

function [values,x]=eigen_values(D,toggle)
global Y b k1 a umax km x3max;
if toggle==1
    x=[0;20;0];
else
    x=steadystate(D);
end
u=(umax*(1-(x(3)/x3max))*x(2))/(km+x(2)+(k1*(x(2)^2)));
dudx2=((km-(k1*(x(2)^2)))*umax*(1-(x(3)/x3max)))/((km+x(2)+(k1*(x(2)^2)))^2);
dudx3=-((umax*x(2))/(x3max*(km+x(2)+(k1*(x(2)^2)))));
r1=[u-D,x(1)*dudx2,x(1)*dudx3];
r2=[-u/Y,-(D+((x(1)/Y)*dudx2)),-(x(1)*dudx3)/Y];
r3=[(a*u)+b,x(1)*((a*dudx2)),(x(1)*((a*dudx3)))-D];
J=[r1;r2;r3];
[~,e_values]=eig(J);
values=[e_values(1,1),e_values(2,2),e_values(3,3)];
end

function dxdt=diff_eq(t,x)
global Y b k1 a umax km x2f Dinit x3max;
u=(umax*(1-(x(3)/x3max))*x(2))/(km+x(2)+(k1*(x(2)^2)));
dx1=(u-Dinit)*x(1);
dx2=(Dinit*(x2f-x(2)))-((u*x(1))/Y);
dx3=(((a*u)+b)*x(1))-(Dinit*x(3));
dxdt=[dx1;dx2;dx3];
end