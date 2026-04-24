%% avr_new
% Compare the generator model without PI and with PI voltage control.
% This script:
% 1) simulates both models from the same plant initial condition,
% 2) computes terminal voltage V(t),
% 3) finds the equilibrium points with fsolve,
% 4) numerically linearizes each model,
% 5) prints the eigenvalues,
% 6) plots root loci of the two linearized plants.
%
% Root-locus convention used here:
% - input  : external signal u added to the Efd equation
% - output : omega (the 4th state)
% This follows the user's requirement: root locus from u to w.

clear; clc; close all;

%% simulation setup
par = get_par();
tspan = [0 200];

% same plant-state initial condition for both cases
x0_nopi = [1.1; 0.9; 0.1; 0.01];
x0_pi   = [1.1; 0.9; 0.1; 0.01; 0.0];
u0 = 0;

opts_ode = odeset('RelTol',1e-9,'AbsTol',1e-9);

%% time-domain simulation
[t_nopi, x_nopi] = ode45(@(t,x) f_nopi_u(t,x,par,u0), tspan, x0_nopi, opts_ode);
[t_pi,   x_pi]   = ode45(@(t,x) f_pi(t,x,par,u0),      tspan, x0_pi,   opts_ode);

V_nopi = calc_V(x_nopi(:,2), x_nopi(:,3), par);
V_pi   = calc_V(x_pi(:,2),   x_pi(:,3),   par);

%% equilibrium points (start from end points of simulation)
opts_fsolve = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-12, ...
    'MaxIterations',1000, ...
    'MaxFunctionEvaluations',5000);

xe_nopi = fsolve(@(x) f_nopi_u(0,x,par,u0), x_nopi(end,:).', opts_fsolve);
xe_pi   = fsolve(@(x) f_pi(0,x,par,u0),      x_pi(end,:).',   opts_fsolve);

Ve_nopi = calc_V(xe_nopi(2), xe_nopi(3), par);
Ve_pi   = calc_V(xe_pi(2),   xe_pi(3),   par); 

%% numerical linearization
h = 1e-6;
A_nopi = numerical_jacobian(@(x) f_nopi_u(0,x,par,u0), xe_nopi, h);
A_pi   = numerical_jacobian(@(x) f_pi(0,x,par,u0),      xe_pi,   h);

B_nopi = numerical_input_jacobian(@(u) f_nopi_u(0,xe_nopi,par,u), u0, length(xe_nopi), h);
B_pi   = numerical_input_jacobian(@(u) f_pi(0,xe_pi,par,u),        u0, length(xe_pi),   h);

% Output 1: terminal voltage V (for reporting only)
C_V_nopi = numerical_output_jacobian(@(x) calc_V(x(2), x(3), par), xe_nopi, h);
C_V_pi   = numerical_output_jacobian(@(x) calc_V(x(2), x(3), par), xe_pi,   h);
D_V_nopi = 0;
D_V_pi   = 0;

% Output 2: omega = 4th state (for root locus)
C_w_nopi = [0 0 0 1];  
C_w_pi   = [0 0 0 1 0];
D_w_nopi = 0;
D_w_pi   = 0;

eig_nopi = eig(A_nopi);
eig_pi   = eig(A_pi);

sysV_nopi = ss(A_nopi, B_nopi, C_V_nopi, D_V_nopi);
sysV_pi   = ss(A_pi,   B_pi,   C_V_pi,   D_V_pi);

sysW_nopi = ss(A_nopi, B_nopi, C_w_nopi, D_w_nopi);
sysW_pi   = ss(A_pi,   B_pi,   C_w_pi,   D_w_pi);

%% display results
format long g;
disp('================ WITHOUT PI ================');
disp('Equilibrium x_e = [Efd; e; delta; omega]');
disp(xe_nopi);
fprintf('Equilibrium V = %.12f\n', Ve_nopi);
disp('A matrix =');
disp(A_nopi);
disp('B matrix (u -> states) =');
disp(B_nopi);
disp('C_V matrix (states -> V) =');
disp(C_V_nopi);
disp('C_w matrix (states -> omega) =');
disp(C_w_nopi);
disp('Eigenvalues =');
disp(eig_nopi);

fprintf('\n');
disp('================= WITH PI =================');
disp('Equilibrium x_e = [Efd; e; delta; omega; W1]');
disp(xe_pi);
fprintf('Equilibrium V = %.12f\n', Ve_pi);
disp('A matrix =');
disp(A_pi);
disp('B matrix (u -> states) =');
disp(B_pi);
disp('C_V matrix (states -> V) =');
disp(C_V_pi);
disp('C_w matrix (states -> omega) =');
disp(C_w_pi);
disp('Eigenvalues =');
disp(eig_pi);

%% plots: compare common states and V
figure('Name','Common states and terminal voltage','Color','w');

subplot(3,2,1);
plot(t_nopi, x_nopi(:,1), 'LineWidth', 1.4); hold on;
plot(t_pi,   x_pi(:,1),   '--', 'LineWidth', 1.4);
grid on; xlabel('Time (s)'); ylabel('E_{fd}');
legend('No PI','With PI','Location','best');
title('Field voltage state');

subplot(3,2,2);
plot(t_nopi, x_nopi(:,2), 'LineWidth', 1.4); hold on;
plot(t_pi,   x_pi(:,2),   '--', 'LineWidth', 1.4);
grid on; xlabel('Time (s)'); ylabel('e');
legend('No PI','With PI','Location','best');
title('Internal emf state');

subplot(3,2,3); 
plot(t_nopi, x_nopi(:,3), 'LineWidth', 1.4); hold on;
plot(t_pi,   x_pi(:,3),   '--', 'LineWidth', 1.4);
grid on; xlabel('Time (s)'); ylabel('\delta (rad)');
legend('No PI','With PI','Location','best');
title('Rotor angle');

subplot(3,2,4);
plot(t_nopi, x_nopi(:,4), 'LineWidth', 1.4); hold on;
plot(t_pi,   x_pi(:,4),   '--', 'LineWidth', 1.4);
grid on; xlabel('Time (s)'); ylabel('\omega');
legend('No PI','With PI','Location','best');
title('Speed deviation');

subplot(3,2,5);
plot(t_nopi, V_nopi, 'LineWidth', 1.4); hold on;
plot(t_pi,   V_pi,   '--', 'LineWidth', 1.4);
yline(par.Vstar, ':', 'V^*', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('V');
legend('No PI','With PI','V^*','Location','best');
title('Terminal voltage');

subplot(3,2,6);
plot(t_pi, x_pi(:,5), '--', 'LineWidth', 1.4);
grid on; xlabel('Time (s)'); ylabel('W_1');
legend('With PI','Location','best');
title('Integrator state');

%% optional: zoom on first 10 seconds
idx_nopi = t_nopi <= 10;
idx_pi   = t_pi   <= 10;
figure('Name','Zoom: first 10 seconds','Color','w');

subplot(2,1,1);
plot(t_nopi(idx_nopi), x_nopi(idx_nopi,4), 'LineWidth', 1.4); hold on;
plot(t_pi(idx_pi),     x_pi(idx_pi,4),   '--', 'LineWidth', 1.4);
grid on; xlabel('Time (s)'); ylabel('\omega');
legend('No PI','With PI','Location','best');
title('Speed deviation (0-10 s)');

subplot(2,1,2);
plot(t_nopi(idx_nopi), V_nopi(idx_nopi), 'LineWidth', 1.4); hold on;
plot(t_pi(idx_pi),     V_pi(idx_pi),   '--', 'LineWidth', 1.4);
yline(par.Vstar, ':', 'V^*', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('V');
legend('No PI','With PI','V^*','Location','best');
title('Terminal voltage (0-10 s)');

%% root locus: u -> omega
figure('Name','Root locus comparison: u to omega','Color','w');
subplot(1,2,1);
rlocus(sysW_nopi);
grid on;
title('Root locus: no PI plant (u -> \omega)');

subplot(1,2,2);
rlocus(sysW_pi);
grid on;
title('Root locus: PI plant (u -> \omega)');

%% transfer functions
fprintf('\nTransfer function from u to omega (no PI):\n');
zpk(sysW_nopi)

fprintf('\nTransfer function from u to omega (with PI):\n');
zpk(sysW_pi)

fprintf('\nTransfer function from u to V (no PI):\n');
zpk(sysV_nopi)

fprintf('\nTransfer function from u to V (with PI):\n');
zpk(sysV_pi)

%% -------- local functions --------
function par = get_par()
    par.Td0p = 0.1;
    par.xdp  = 0.3;
    par.xd   = 1.8;
    par.xT   = 0.1;
    par.x    = par.xd + par.xdp + par.xT;
    par.xL   = 1.9;
    par.alpha = par.xL/(par.x + par.xL);
    par.beta  = par.x /(par.x + par.xL);
    par.K1 = 1.5;
    par.K2 = 6;
    par.K3 = 0.5;
    par.Td = 0.05;
    par.a = par.xd/par.xdp;
    par.b = (par.xd-par.xdp)/par.xdp;
    par.Vstar = 0.9;
    par.M  = 0.2;
    par.Pm = 0.1;
    par.d  = 0.28;
end

function dx = f_nopi_u(~,x,par,u)
    efd = x(1);
    e   = x(2);
    del = x(3);
    om  = x(4);

    [V, theta] = V_theta(e, del, par);

    dx = zeros(4,1);
    dx(1) = (-par.K1*efd - par.K2*(V - par.Vstar) + u)/par.Td0p;
    dx(2) = (-par.a*e + par.b*V*cos(del-theta) + efd)/par.Td;
    dx(3) = om;
    dx(4) = (par.Pm - par.d*om - e*V*sin(del-theta)/par.x)/par.M;
end

function dx = f_pi(~,x,par,u)
    efd = x(1);
    e   = x(2);
    del = x(3);
    om  = x(4);
    W1  = x(5);

    [V, theta] = V_theta(e, del, par);

    dx = zeros(5,1);
    dx(1) = (-par.K1*efd - par.K2*(V - par.Vstar) - par.K3*W1 + u)/par.Td0p;
    dx(2) = (-par.a*e + par.b*V*cos(del-theta) + efd)/par.Td;
    dx(3) = om;
    dx(4) = (par.Pm - par.d*om - e*V*sin(del-theta)/par.x)/par.M;
    dx(5) = V - par.Vstar;
end

function [V, theta] = V_theta(e, del, par)
    V = sqrt(par.alpha^2*e.^2 + 2*par.alpha*par.beta*e.*cos(del) + par.beta^2);
    theta = atan((par.alpha*e.*sin(del))./(par.alpha*e.*cos(del) + par.beta));
end

function V = calc_V(e, del, par)
    V = sqrt(par.alpha^2*e.^2 + 2*par.alpha*par.beta*e.*cos(del) + par.beta^2);
end

function A = numerical_jacobian(fun, xeq, h)
    n = length(xeq);
    A = zeros(n,n);
    for i = 1:n
        dx = zeros(n,1);
        dx(i) = h;
        fp = fun(xeq + dx);
        fm = fun(xeq - dx);
        A(:,i) = (fp - fm)/(2*h);
    end
end

function B = numerical_input_jacobian(fun_u, ueq, n, h)
    fp = fun_u(ueq + h);
    fm = fun_u(ueq - h);
    B = (fp - fm)/(2*h);
    B = reshape(B, [n, 1]);
end

function C = numerical_output_jacobian(fun_y, xeq, h)
    n = length(xeq);
    C = zeros(1,n);
    for i = 1:n
        dx = zeros(n,1);
        dx(i) = h;
        yp = fun_y(xeq + dx);
        ym = fun_y(xeq - dx);
        C(i) = (yp - ym)/(2*h);
    end
end
