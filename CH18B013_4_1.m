%Parameters
global D L H ff v0 A g rho c deltaz n
D = 0.5;
L = 600;
H = 100;
ff = 5e-3;
g = 9.81;
v0 = ((g*H)/(0.5+(2*ff*L/D)))^0.5;
A = (pi*(D^2))/4;
rho = 1000;
c = 1200;
deltaz = 100;
n = (L/deltaz)+1;
tfinal = 0;
z = 0:deltaz:L;

%Initial conditions and ode solver
P0 = rho*g*H*ones(1,n);
m0 = zeros(1,n);
y0 = [P0(2:n),m0(1:n-1)];
[t,y] = ode23s(@pipedynfun,[10:-0.01:tfinal],y0);
P1 = rho*g*H*ones(length(t),1);
mn = rho*A*v0*(1-(t/10));

%Plot of Pressure vs time
figure
plot(t,P1)
hold on
for i = 1:n-1
    plot(t,y(:,i))
    if i == n-1
        hold off
    end
    hold on
end
title('Plot of Pressure vs time')
xlabel('Time (s)')
ylabel('Pressure (Pa)')
legend('P1','P2','P3','P4','P5','P6','P7')

%Plot of velocity vs time
figure
for i = n:2*(n-1)
    plot(t,y(:,i)/(rho*A))
    hold on
end
plot(t,mn/(rho*A))
title('Plot of velocity vs time')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
legend('v1','v2','v3','v4','v5','v6','v7')

%Explanation:
%The transient profiles for Pressure and velocity have sinusoidal nature
%similar to the the analysis done when the height of the water level in
%tank is varied sinusoidally. We assume that at the end of 10 seconds, the
%water in the pipe stops flowing completely and has a static pressure
%distribution. 

%Function definition
function f = pipedynfun(t,y)
global D H ff v0 A g rho c deltaz n
P = y(1:n-1);
m = y(n:2*n-2);
f = zeros(2*n-2,1);
P1 = rho*g*H;
mn = rho*A*v0*(1-(t/10));
P = [P1;P];
m = [m;mn];
for i = 2:n-1
    f(i-1,1) = (c^2/A)*(m(i-1)-m(i+1))/(2*deltaz);
    f(i+n-1,1) = A*(P(i-1)-P(i+1))/(2*deltaz)-A*(2*ff/(A^2*D))*m(i)*abs(m(i))/rho;
end
f(n-1,1) = (c^2/A)*(m(n-1)-m(n))/(deltaz);
f(n,1) = A*(P(1)-P(2))/(deltaz)-A*(2*ff/(A^2*D))*m(1)*abs(m(1))/rho;
end