global R Aw E0 A2 E2 A34 E34 A3g E3g Mwmax Mwmin mw0 A B C D E F m2vec Tvec tvec rhow delH r
R=8.314;
Aw=103500;
E0=87500;
A=0.10;B=0.41;C=0.22;D=0.27;E=0.5;F=0.5;
A2=1.18;
E2=45000;
A34=2e11;
E34=251000;
A3g=3e13;
E3g=301500;
Mwmax=0.45;Mwmin=0.17;
rhow=510;
mw0=0.02326412192;%For beech wood, mw0=(4.549458131e-3)*(3872/529).For pine wood, mw0=0.02326412192 
delH=-125000;
r=0.011;
m2vec=[];Tvec=[];tvec=[];

%Beech, (0.05K/s,15000s,275K), (1K/s,700s,290K)
%Pine, (0.05K/s,15000s,275K), (1K/s,700s,290K)
[time,mass]=ode45(@diff_eq,[0:50:200],[275;275;275;mw0;0;0;0]);
char=mass(:,7);
wood=mass(:,4);
Temperature_centre=mass(:,1);
Temperature_surface=mass(:,3);
gasflow=((E*A2*mass(:,5)).*(exp((-E2/R)*(1./mass(:,1)))))+((A3g*mass(:,6)).*(exp((-E3g/R)*(1./mass(:,1)))))+((A*Aw*mass(:,4)).*(exp((-E0/R)*(1./mass(:,1)))));
gas=cumtrapz(time,gasflow);
tar=(mw0*ones(1,length(time)))-(wood+gas+char);

figure
%plot(time,Temperature)
%hold on
plot(time,Temperature_centre)
hold on
plot(time,Temperature_surface)
hold off
legend('T in','T ext')

function y=fun(m,T)
global R A2 E2
y=m*A2*exp((-E2)/(R*T));
end

function dmdt=diff_eq(t,m)
global R Aw E0 A2 E2 A34 E34 A3g E3g Mwmax Mwmin mw0 A C D E F m2vec Tvec tvec rhow delH r
tvec=[tvec,t];
Tvec=[Tvec,m(1)];
m2vec=[m2vec,m(3)];
funvec=zeros(1,length(t));
for i=1:length(tvec)
    funvec(i)=fun(m2vec(i),Tvec(i));
end
mdotg=((E*A2*m(5))*(exp((-E2/R)*(1/m(1)))))+((A3g*m(6))*(exp((-E3g/R)*(1/m(1)))))+((A*Aw*m(4))*(exp((-E0/R)*(1/m(1)))));
Rp=m(4)*Aw*exp((-E0)/(R*m(1)));
mw=m(4)/(m(4)+m(5)+m(6)+m(7));
mc=m(7)/(m(4)+m(5)+m(6)+m(7));
k=(mc*0.105)+(mw*((rhow*0.1941)+0.0186));
Cp1=(mc*((0.36*m(1))+1390))+(mw*(1112.0+(4.85*(m(1)-273.15))));
Cp2=(mc*((0.36*m(2))+1390))+(mw*(1112.0+(4.85*(m(2)-273.15))));
Cp3=(mc*((0.36*m(3))+1390))+(mw*(1112.0+(4.85*(m(3)-273.15))));
d2Tdr2=((2*m(3))-(4*m(2))+(2*m(1)))/(r^2);
dTdr1=(m(2)-m(1))/(r/2);
dTdr2=(m(3)-m(1))/(r);
dTdr3=(m(3)-m(2))/(r/2);
dmdt1=((k/(rhow*Cp1))*(d2Tdr2+(dTdr1/1e-5)))-((mdotg/(2*pi*rhow))*dTdr1)+((Rp*delH)/(rhow*Cp1));
dmdt2=((k/(rhow*Cp2))*(d2Tdr2+(dTdr2/(r/2))))-((mdotg/(2*pi*rhow))*dTdr2)+((Rp*delH)/(rhow*Cp2));
dmdt3=((k/(rhow*Cp3))*(d2Tdr2+(dTdr3/r)))-((mdotg/(2*pi*rhow))*dTdr3)+((Rp*delH)/(rhow*Cp3));
dmdt4=-m(4)*Aw*exp((-E0)/(R*m(1)));
f2=1-(trapz(tvec,funvec,2)/(C*mw0));
m20=((mw0*C)/((Mwmax-Mwmin)*1e3))*(((log(87060)/299)^(1/0.59))/0.59)*(m(1)^((1/0.59)-1))*dmdt1;
dmdt5=-(m(5)*A2*exp((-E2)/(R*m(1))))+(f2*m20);
dmdt6=(D*m(4)*Aw*exp((-E0)/(R*m(1))))+(F*m(5)*A2*exp((-E2)/(R*m(1))))-(m(6)*A34*exp((-E34)/(R*m(1))))-(m(6)*A3g*exp((-E3g)/(R*m(1))));
dmdt7=m(6)*A34*exp((-E34)/(R*m(1)));
dmdt=[dmdt1;dmdt2;dmdt3;dmdt4;dmdt5;dmdt6;dmdt7];
end