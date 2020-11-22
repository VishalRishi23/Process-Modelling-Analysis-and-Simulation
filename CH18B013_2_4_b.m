%ode15s with combined approach using mass matrix
global pL
pL=1000; %Density of fluid (kg/m3)

%Solver for time less than 416s at which P<P0
syms VL(t) mdv(t) mv(t) VV(t) TL(t) P(t) V Qin Qh Tin rhoL CpL Tr H MW R a b A B C;
eq1 = diff(VL,1)==Qin-(mdv/rhoL);
eq2 = diff(mv,1)==mdv;
eq3 = mdv==((CpL*(Tin-Tr)*rhoL*Qin)+Qh)/((CpL*(TL-Tr))+H);
eq4 = (((P+((a*(mv^2))/((MW*VV)^2)))*(VV-((mv*b)/MW)))-(mv*R*TL/MW))==0;
eq5 = ((133.322*(10^(A-(B/(TL+C-273.15)))))-P)==0;
eq6 = VV==V-VL;
eqns = [eq1 eq2 eq3 eq4 eq5 eq6]; %Equations
vars=[VL(t);mdv(t);mv(t);VV(t);TL(t);P(t)]; %Variables
[DAEs,DAEvars]=reduceDAEIndex(eqns,vars);
[M,f]=massMatrixForm(DAEs,DAEvars); %Mass matrix
M=odeFunction(M,DAEvars);
f=odeFunction(f,DAEvars,V,Qin,Qh,Tin,rhoL,CpL,Tr,H,MW,R,a,b,A,B,C);
V=5.663;Qin=2.816E-4;Qh=5861;Tin=294.26;rhoL=1000;CpL=4186;Tr=273.15;
H=1984480;MW=0.018;R=8.314;a=0.5537;b=0.03049E-3;A=8.14019;B=1810.94;C=244.485; %Parameters
F=@(t, Y) f(t,Y,V,Qin,Qh,Tin,rhoL,CpL,Tr,H,MW,R,a,b,A,B,C);
y0=[2.832;0.0116;9.072;2.831;434.0861;628462.568318097];
yp0=zeros(6,1);
opt=odeset('Mass', M, 'InitialSlope', yp0,'RelTol', 10.0^(-7), 'AbsTol' , 10.0^(-7));
[t1,y1] = ode15s(F,[0;416],y0,opt);

%Solver for time greater than 416s at which P>=P0
syms VL(t) mdv(t) mv(t) mdvout(t) VV(t) TL(t) P(t) V Qin Qh Tin P0 rhoL CpL Tr H Kv MW R a b A B C;
eq1 = diff(VL,1)==Qin-(mdv/rhoL);
eq2 = diff(mv,1)==mdv-mdvout;
eq3 = mdv==((CpL*(Tin-Tr)*rhoL*Qin)+Qh)/((CpL*(TL-Tr))+H);
eq4 = (((P+((a*(mv^2))/((MW*VV)^2)))*(VV-((mv*b)/MW)))-(mv*R*TL/MW))==0;
eq5 = ((133.322*(10^(A-(B/(TL+C-273.15)))))-P)==0;
eq6 = VV==V-VL;
eq7 = mdvout==Kv*(((P-P0)*P)^(0.5));
eqns = [eq1 eq2 eq3 eq4 eq5 eq6 eq7]; %Equations
vars=[VL(t);mdv(t);mv(t);mdvout(t);VV(t);TL(t);P(t)]; %Variables
[DAEs,DAEvars]=reduceDAEIndex(eqns,vars);
[M,f]=massMatrixForm(DAEs,DAEvars); %Mass matrix
M=odeFunction(M,DAEvars);
f=odeFunction(f,DAEvars,V,Qin,Qh,Tin,P0,rhoL,CpL,Tr,H,Kv,MW,R,a,b,A,B,C);
V=5.663;Qin=2.816E-4;Qh=5861;Tin=294.26;P0=1034200;rhoL=1000;CpL=4186;Tr=273.15;
H=1984480;Kv=1.101E-7;MW=0.018;R=8.314;a=0.5537;b=0.03049E-3;A=8.14019;B=1810.94;C=244.485; %Parameters
F=@(t, Y) f(t,Y,V,Qin,Qh,Tin,P0,rhoL,CpL,Tr,H,Kv,MW,R,a,b,A,B,C);
y0=[2.9461;0.0112;13.8434;0.0036;2.7169;454.7631;[1035230.38149435]];
yp0=zeros(7,1);
opt=odeset('Mass', M, 'InitialSlope', yp0,'RelTol', 10.0^(-7), 'AbsTol' , 10.0^(-7));
[t2,y2] = ode15s(F,[422.3610;1500],y0,opt);

%Plots
figure
plot(t1,y1(:,5),t2,y2(:,6));title('Temperature in the boiler');xlabel('time');ylabel('T (K)');
figure
plot(t1,y1(:,6),t2,y2(:,7));title('Pressure in the boiler');xlabel('time');ylabel('P (Pa)');
figure
plot(t1,y1(:,2),t2,y2(:,2));title('Vapour mass flow rate leaving the liquid');xlabel('time');ylabel('mdv (kg/s)');
figure
plot(t2,y2(:,4));title('Vapour mass flow rate leaving the boiler');xlabel('time');ylabel('mdvout (kg/s)');
figure
plot(t1,y1(:,1),t2,y2(:,1));title('Volume of liquid in the boiler');xlabel('time');ylabel('VL (m3)');
figure
plot(t1,y1(:,4),t2,y2(:,5));title('Volume of vapour in the boiler');xlabel('time');ylabel('VV (m3)');
figure
plot(t1,pL*y1(:,1),t2,pL*y2(:,1));title('Mass of liquid in the boiler');xlabel('time');ylabel('ml (kg)');
figure
plot(t1,y1(:,3),t2,y2(:,3));title('Mass of vapour in the boiler');xlabel('time');ylabel('mv (kg)');