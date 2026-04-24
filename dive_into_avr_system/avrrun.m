close all; clc; clear;

tspan=[0 10];
y0=[1.1; 0.9; 0.1; 0.01; 0;0 ]; 


[tout, y] = ode45(@(t,y) avr(t,y), tspan, y0);



figure;

plot(tout,y(:,4)); grid on; ylabel('\omega'); xlabel('Time (s)');


tspan = [0 10];

% same initial condition for the plant states
y0_nopss = [1.1; 0.9; 0.1; 0.01];
y0_pss   = [1.1; 0.9; 0.1; 0.01; 0; 0];

% without PSS
[t1, y1] = ode45(@(t,y) avr_nopss(t,y), tspan, y0_nopss);

% with PSS
[t2, y2] = ode45(@(t,y) avr(t,y), tspan, y0_pss);

% plot omega on the same figure
figure;
plot(t1, y1(:,4), 'LineWidth', 1.5); hold on;
plot(t2, y2(:,4), '--', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('\omega');
legend('without PSS','with PSS');
title('\omega response');

% % voltage
% e   = y(:,2);
% del = y(:,3);
% xdp=0.3;
% xd=1.8;
% xT=0.001;
% x=xd+xdp+xT;
% xL=1.9;
% alpha=xL/(x+xL);
% beta = x/(x+xL);
% 
% V = sqrt(alpha^2 .* e.^2 + 2*alpha*beta .* e .* cos(del) + beta^2);
% figure;
% plot(tout,V,'LineWidth',1.2);     grid on; ylabel('V'); xlabel('Time (s)');
% 
% 
% %% linearize
% y0 = [1.1; 0.9; 0.1; 0];    % omega=0
% 
% % % equili point 
% % opts = optimoptions('fsolve', ...
% %     'Display','iter', ...
% %     'FunctionTolerance',1e-12, ...
% %     'StepTolerance',1e-12);
% % 
% % y0 = fsolve(@(y) avr(0,y,u0), y0_guess, opts);
% % 
% % disp(y0);
% 
% % linearization
% n = length(y0);
% A = zeros(n,n);
% B = zeros(n,1);
% 
% h = 1e-6;
% 
% % A
% for i = 1:n
%     dy = zeros(n,1);
%     dy(i) = h;
%     f_plus  = avr(0, y0 + dy, u0);
%     f_minus = avr(0, y0 - dy, u0);
%     A(:,i) = (f_plus - f_minus)/(2*h);
% end
% 
% % B 
% f_plus  = avr(0, y0, u0 + h);
% f_minus = avr(0, y0, u0 - h);
% B = (f_plus - f_minus)/(2*h);
% 
% C = [0 0 0 1];
% D = 0;
% 
% % ss model
% sys = ss(A,B,C,D);
% 
% 
% A
% B
% C
% D
% 
% %% ss--tf
% G = tf(sys)
% 
% [num, den] = tfdata(G, 'v');
% 
% pole(G)
% zero(G)


%% rl
s = tf('s');
Yi = 1/s;
% rltool(G); % we want to eliminate the zeros in G(s), so C(s) has 1/s term

%% C(s)  tf to ss
Cc = 61.96*(1 + 1.3*s + (0.98*s)^2)/(s*(1 + s/179.6));




Cc

%% Convert controller Cnew(s) to state-space
sysC = ss(Cc);

[A1,B1,C1,D1] = ssdata(sysC);

disp('A1 = ')
disp(A1)

disp('B1 = ')
disp(B1)

disp('C1 = ')
disp(C1)

disp('D1 = ')
disp(D1)