global R Aw E0 A2 E2 A34 E34 A3g E3g Mwmax Mwmin mw0 A B C D E F m2vec Tvec tvec dTdt
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
mw0=(4.549458131e-3)*(3872/529); %For beech wood, mw0=(4.549458131e-3)*(3872/529).For pine wood, mw0=0.02326412192 
dTdt=1.0;
m2vec=[];Tvec=[];tvec=[];

%Beech, (0.05K/s,15000s,275K), (1K/s,700s,290K)
%Pine, (0.05K/s,15000s,275K), (1K/s,700s,290K)
[time,mass]=ode45(@diff_eq,[0:50:700],[290;mw0;0;0;0]);
char=mass(:,5);
wood=mass(:,2);
Temperature=mass(:,1);
gasflow=((E*A2*mass(:,3)).*(exp((-E2/R)*(1./mass(:,1)))))+((A3g*mass(:,4)).*(exp((-E3g/R)*(1./mass(:,1)))))+((A*Aw*mass(:,2)).*(exp((-E0/R)*(1./mass(:,1)))));
gas=cumtrapz(time,gasflow);
tar=(mw0*ones(1,length(time)))-(wood+gas+char);

figure
%plot(time,Temperature)
%hold on
plot(time,gas/mw0)
hold on
plot(time,tar/mw0)
hold off
legend('gas','tar')

function y=fun(m,T)
global R A2 E2
y=m*A2*exp((-E2)/(R*T));
end

function dmdt=diff_eq(t,m)
global R Aw E0 A2 E2 A34 E34 A3g E3g Mwmax Mwmin mw0 C D F m2vec Tvec tvec dTdt
tvec=[tvec,t];
Tvec=[Tvec,m(1)];
m2vec=[m2vec,m(3)];
funvec=zeros(1,length(t));
for i=1:length(tvec)
    funvec(i)=fun(m2vec(i),Tvec(i));
end
dmdt1=dTdt;
dmdt2=-m(2)*Aw*exp((-E0)/(R*m(1)));
f2=1-(trapz(tvec,funvec,2)/(C*mw0));
m20=((mw0*C)/((Mwmax-Mwmin)*1e3))*(((log(87060)/299)^(1/0.59))/0.59)*(m(1)^((1/0.59)-1))*dTdt;
dmdt3=-(m(3)*A2*exp((-E2)/(R*m(1))))+(f2*m20);
dmdt4=(D*m(2)*Aw*exp((-E0)/(R*m(1))))+(F*m(3)*A2*exp((-E2)/(R*m(1))))-(m(4)*A34*exp((-E34)/(R*m(1))))-(m(4)*A3g*exp((-E3g)/(R*m(1))));
dmdt5=m(4)*A34*exp((-E34)/(R*m(1)));
dmdt=[dmdt1;dmdt2;dmdt3;dmdt4;dmdt5];
end