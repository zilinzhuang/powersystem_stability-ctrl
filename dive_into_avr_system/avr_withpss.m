function yd=avr(t,y)

yd=zeros(6,1);

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
Td=5;
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
gamma1=y(5);
gamma2=y(6);



V=sqrt(alpha^2*e^2+2*alpha*beta*e*cos(del)+beta^2);
theta=atan(alpha*e*sin(del)/(alpha*e*cos(del)+beta));
u=-(-1860*gamma1 +10.9*gamma2)/(1+1.0687e+04);

yd(1)=1/Td0p*(-y(1)*K1-K2*(V-Vstar)+u);
yd(2)=1/Td*(-a*y(2)+b*V*cos(del-theta)+y(1));
yd(3)=y(4);
yd(4)=1/M*(Pm-d*y(4)-e*V*sin(del-theta)/x);
yd(5)= -179.6*y(5)-512*2*y(4);
yd(6)= y(5);


