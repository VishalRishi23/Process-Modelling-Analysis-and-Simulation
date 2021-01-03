%Parameters
global D L a v0 A g Ti deltaz n
D = 0.1;
L = ;
a = 10e-4;
g = 9.81;
v0 = 0.5;
A = (pi*(D^2))/4;
Ti = 300;
deltaz = 1/3;
n = (L/deltaz)+1;
z = 0:deltaz:L;

%Initial conditions and bvp solver
rmesh = linspace(0,0.05,50);
solinit = bvpinit(rmesh,zeros(1,2*n-2));
sol = bvp4c(@bvpfun, @bcfcn, solinit);
T1 = Ti*ones(length(sol.x),1);

%Plot of Temperature vs radial distance
figure
plot(sol.x,T1)
hold on
for i = 1:n-1
    plot(sol.x,sol.y(i,:))
    if i == n-1
        hold off
    end
    hold on
end
title('Plot of Temperature vs radial distance')
xlabel('Radial distance (m)')
ylabel('Temperature (K)')
legend('T1','T2','T3','T4','T5','T6','T7')


%Function definitions
function f = bvpfun(r,y)
global D a v0 Ti deltaz n
v = v0*(1-(((r/D)^2)*4));
T = y(1:n-1);
Td = y(n:2*n-2);
f = zeros(2*n-2,1);
T1 = Ti;
Td1 = 0;
T = [T1;T];
Td = [Td1;Td];
for i = 2:n
    f(i-1,1) = Td(i);
    f(i+n-2,1) = ((v/a)*((T(i)-T(i-1))/deltaz))-(Td(i)/(r+10e-20));
end
end

function res = bcfcn(ya,yb)
%res = [ya(7);yb(7)-10;ya(8);yb(8)-10;ya(9);yb(9)-10;ya(10);yb(10)-10;ya(11);yb(11)-10;ya(12);yb(12)-10];
res = [ya(1)-300;ya(2)-300;ya(3)-300;ya(4)-300;ya(5)-300;ya(6)-300;yb(7)-10;yb(8)-10;yb(9)-10;yb(10)-10;yb(11)-10;yb(12)-10];
end