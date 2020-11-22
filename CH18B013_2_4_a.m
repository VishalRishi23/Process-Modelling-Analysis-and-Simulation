%ode15s with nested approach
global V Qin Qh Tin P0 pL CpL Tr H Kv MW R a b A B C;
V=5.663; %Total Volume (m3)
Qin=2.816E-4; %Inlet flow rate (m3/s)
Qh=5861; %Heat transfer rate (J/s)
Tin=294.26; %Inlet temperature (K)
P0=1034200; %Downstream pressure (Pa)
pL=1000; %Density of fluid (kg/m3)
CpL=4186; %Specific heat capacity (J/K/kg)
Tr=273.15; %Reference temperature (K)
H=1984480; %Latent heat (J/kg)
Kv=1.101E-7;
MW=0.018; %Molecular weight (kg/mol)
R=8.314;  %Universal gas constant (J/K/mol)
a=0.5537;
b=0.03049E-3;
A=8.14019;
B=1810.94;
C=244.485;
[t,v]=ode15s(@diff_eq,[0;1500],[2.832;9.072]);
VLvec=v(:,1); %Volume of liquid (m3)
mvvec=v(:,2); %Mass of vapour (kg)
[~,VV,P,TL,mdvout,mdv]=cellfun(@(t,v) diff_eq(t,v.'),num2cell(t),num2cell(v,2),'uni',0);
VVvec=cell2mat(VV); %Volume of vapour (m3)
Pvec=cell2mat(P); %Pressure in the boiler (Pa) 
TLvec=cell2mat(TL); %Temperature in the boiler (K)
mdvoutvec=cell2mat(mdvout); %Vapour flow rate leaving the boiler (kg/s)
mdvvec=cell2mat(mdv); %Vapour flow rate leaving the liquid (kg/s)
mlvec=pL*VLvec; %Mass of liquid (kg)

%Plots
figure
plot(t,TLvec);title('Temperature in the boiler');xlabel('time');ylabel('T (K)');
figure
plot(t,Pvec);title('Pressure in the boiler');xlabel('time');ylabel('P (Pa)');
figure
plot(t,mdvvec);title('Vapour mass flow rate leaving the liquid');xlabel('time');ylabel('mdv (kg/s)');
figure
plot(t,mdvoutvec);title('Vapour mass flow rate leaving the boiler');xlabel('time');ylabel('mdvout (kg/s)');
figure
plot(t,VLvec);title('Volume of liquid in the boiler');xlabel('time');ylabel('VL (m3)');
figure
plot(t,VVvec);title('Volume of vapour in the boiler');xlabel('time');ylabel('VV (m3)');
figure
plot(t,mlvec);title('Mass of liquid in the boiler');xlabel('time');ylabel('ml (kg)');
figure
plot(t,mvvec);title('Mass of vapour in the boiler');xlabel('time');ylabel('mv (kg)');

%Function definitions
function [dvdt,VV,P,TL,mdvout,mdv]=diff_eq(t,v)
    global V Qin Qh Tin P0 pL CpL Tr H Kv MW R a b A B C;
    VV=V-v(1);
    options=optimoptions('fsolve');
    prop=fsolve(@(x) solver(x,v(2),VV),[2E+6;500],options);
    P=prop(1);
    TL=prop(2);
    if P>P0
        mdvout=Kv*(((P-P0)*P)^(0.5));
    else
        mdvout=0;
    end
    mdv=((CpL*(Tin-Tr)*pL*Qin)+Qh)/((CpL*(TL-Tr))+H);
    dv1=Qin-(mdv/pL);
    dv2=mdv-mdvout;
    dvdt=[dv1;dv2];
    function f=solver(x,mv,VV)
        f1=(((x(1)+((a*(mv^2))/((MW*VV)^2)))*(VV-((mv*b)/MW)))-(mv*R*x(2)/MW))/1E+6;
        f2=((133.322*(10^(A-(B/(x(2)+C-273.15)))))-x(1))/1E+6;
        f=[f1;f2];
    end
end