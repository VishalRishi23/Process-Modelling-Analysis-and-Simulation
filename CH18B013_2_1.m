global A F0 beta;
A=0.5; %Area of cross section (m^2)
F0=3; %Flow rate (m^3/min)
h0=5; %Initial height (m)
beta=[0.75,1]; %Proportionality constants (m^2.5/min)
[t,h]=ode45(@diff_eq,[0;100],[h0;h0]);

%Plot of h1 and h2 as a function of time
figure
hold on   
plot(t,h(:,1),'bo')
hold on
plot(t,h(:,2),'rs')
title('Plot of h1 and h2 as a function of time')
xlabel('Time t')
ylabel('Height h')
legend('h1','h2')
hold off

%Function definition to solve the ode's
function dhdt=diff_eq(t,h)
    global A F0 beta;
    dh1=(F0/A)-((beta(1)/A)*((h(1)-h(2))^(0.5)));
    dh2=((beta(1)*((h(1)-h(2))^(0.5)))-(beta(2)*(h(2)^(0.5))))/A;
    dhdt=[dh1;dh2];
end

%Explanation:
%The water level in both tanks reaches steady state. The water level in
%tank 2 attains a minimum before reaching the steady state. This is because
%the initial difference in the water level of the two tanks is zero. There
%is no outlet flow from tank 1. Hence, initially, for a short interval, h1
%increases and h2 decreases (because of outlet flow from tank 2). As the
%height difference increases, there is considerable amount of outlet flow
%from tank 1. Eventually, steady state is acieved.