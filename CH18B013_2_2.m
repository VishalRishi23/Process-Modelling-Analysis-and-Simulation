global pCp V Vj Fj Ti UA Tji;
pCp=4111160;
V=0.283168; %Volume of tank (m3)
Vj=0.0283168; %Volume of jaket (m3)
Fj=7E-4; %Jacket flow rate (m3/s)
Ti=283.15; %Inlet tank temperature (K)
UA=5820; %Overall heat transfer coefficient (J/K/s)
Tji=366.483; %Inlet jacket temperature (K)
[t,T]=ode23s(@diff_eq,[0;1500],[Ti;Tji]);

%Plot of T and Tj as a function of time
figure
hold on
plot(t,T(:,1),'b-','MarkerSize',5)
hold on
plot(t,T(:,2),'r-','MarkerSize',5)
title('Plot of T and Tj as a function of time')
xlabel('Time t')
ylabel('Temperature T')
legend('T','Tj')
hold off

%Function definitions
function F=flowrate(t)
    F=((0.2*sin(5*t))+5)*(10^(-4));
end

function dTdt=diff_eq(t,T)
    global pCp V Vj Fj Ti UA Tji;
    dT1=((Ti-T(1))*(flowrate(t)/V))+((T(2)-T(1))*(UA/(pCp*V)));
    dT2=((Tji-T(2))*(Fj/Vj))-((T(2)-T(1))*(UA/(pCp*Vj)));
    dTdt=[dT1;dT2];
end

%Explanation:
%The temperatures of the tank and the jacket reach a steady state.
%The temperature of the tank increases steadily before reaching the steady steady state.
%The jacket temperature attains a minimum before reaching the steady state.
%Due to the initial temperature difference, there is a rapid fall in Tj.
%This can be explained by viewing plot of heat transfer rate vs time. 
%The heat transfer rate falls steeply and suddenly becomes flat. 
%This steep region causes the rapid fall of the jacket temperature. 
%After the steep region, the heat transfer rate becomes regular which helps
%the temperatures in the tank and jacket reach their steady values.