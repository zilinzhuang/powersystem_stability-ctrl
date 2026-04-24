%% avr_withpi

function yd=avr_withpi(t,y,u)

yd=zeros(4,1);

Td0p=0.1;
xdp=0.3;
xd=1.8;
xT=0.001;
x=xd+xdp+xT;
xL=1.9;
alpha=xL/(x+xL);
beta = x/(x+xL);
K1=1.5;
K2=6;
K3=0.5;

Td=0.05;
a=xd/xdp;
b=(xd-xdp)/xdp;
Vstar=0.9;

M=0.2;
Pm=0.1;
d=0.28;

efd=y(1);
e=y(2);
del=y(3);
om=y(4);
W1=y(5);


V=sqrt(alpha^2*e^2+2*alpha*beta*e*cos(del)+beta^2);
theta=atan(alpha*e*sin(del)/(alpha*e*cos(del)+beta));


yd(1)=1/Td0p*(-y(1)*K1-K2*(V-Vstar)-K3*y(5)+u);
yd(2)=1/Td*(-a*y(2)+b*V*cos(del-theta)+y(1));
yd(3)=y(4);
yd(4)=1/M*(Pm-d*y(4)-e*V*sin(del-theta)/x);
yd(5) = V-Vstar;

