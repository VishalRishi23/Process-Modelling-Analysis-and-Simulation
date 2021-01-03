% Parameters
global v rhocp E deltaEr deltaH cAin Tin deltaz zfinal PeH PeM
zfinal = 2; 
tfinal = 6; 
deltaz = 0.01;
n = zfinal/deltaz+1;
cAin = 1;
Tin = 450;
E = 60000;
deltaEr = -10000; 
deltaH = -100000;
rhocp = 800;
v = 0.4;
PeH = 1; 
PeM = 2.5;

% Steady state solution
zmesh = linspace(0,2);
solinit = bvpinit(zmesh,@guess);
sol = bvp4c(@nonisopfrssfun, @bcfun, solinit);
cAss = sol.y(3,:); Tss = sol.y(1,:);
z = sol.x;
figure
subplot(2,1,1), plot(z,cAss,'LineWidth',1.5), xlabel('Distance'); ylabel('Concentration');title('Concentration vs Distance');
subplot(2,1,2), plot(z,Tss,'LineWidth',1.5), xlabel('Distance'); ylabel('Temperature');title('Temperature vs Distance');

% Functions
function f = nonisopfrssfun(z,y)
    global v rhocp E deltaEr deltaH cAin Tin deltaz zfinal PeH PeM
    T = y(1); dTdz = y(2); cA = y(3); dcAdz = y(4);
    k0 = 2; K0 = 1; R = 8.314; T1 = 450;
    k = k0*exp(-(E/R)*(1/T-1/T1));
    kr = K0*exp(-(deltaEr/R)*(1/T-1/T1));
    rA = -1*k*cA/(kr+cA);
    f(1) = dTdz;
    f(2) = (dTdz - deltaH*rA/(rhocp*v))*PeH/zfinal;
    f(3) = dcAdz;
    f(4) = (dcAdz - (rA/v))*PeM/zfinal;
    f = f';
end

function res = bcfun(in, out)
    global v rhocp E deltaEr deltaH cAin Tin deltaz zfinal PeH PeM
    T_in = in(1); dTdz_in = in(2);
    cA_in = in(3); dcAdz_in = in(4);
    T_out = out(1); dTdz_out = out(2);
    cA_out = out(3); dcAdz_out = out(4);
    res = [T_in - Tin; dTdz_out; cA_in - zfinal*dcAdz_in/PeM - cAin; dcAdz_out];
end

function g = guess(z)
    g = [450+(10*z); 10; 1-(0.1*z); -0.1];
end