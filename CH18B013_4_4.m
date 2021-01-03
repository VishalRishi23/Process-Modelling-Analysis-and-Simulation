% Parameters
zinit = 0;
zfin = 2;
tinit = 0;
tfin = 10;
deltaz = 0.001;
deltat = 0.001;
c0 = 1;
v = 0.25;

m = (zfin - zinit)/deltaz + 1;
n = (tfin - tinit)/deltat + 1;
p = v*deltat/deltaz;
zvec = zinit : deltaz : zfin;
tvec = tinit : deltat : tfin;

% Physical conditions
% At t=0, the concentration in the PFR is 0
c(1:m,1) = 0;
% At t=deltat, the tracer liquid is injected at z=0
c(1,2) = c0;

% Calculations using FTBS scheme
for j = 1 : n-1
    for i = 2 : m
        c(i,j+1) = c(i,j) - p*(c(i,j)-c(i-1,j));
    end
end

% Ploting the results
figure
plot(zvec(1:m),c(1:m,1:500:n),'LineWidth',1.5);
title('Pulse flow in PFR');
xlabel('Distance (m)'); 
ylabel('Concentration (mol/L)')

%Explanation:
%The pulse moves forward along the length of the PFR. We obtain peaks and
%it takes some time to inject the liquid. This is why we obtain a peak
%instead of a vertical line. The area under the peak is a measure of the
%amount of tracer liquid. Once injected, the water carries out the liquid
%and the pulse propagates. This is the RTD experiment.