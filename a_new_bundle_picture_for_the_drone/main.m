%% Simulation.
global chart_pre
% Desired state.
Xd_fun = @(t) [0 + t - t; % x
               sin(t);
               1 - 0.5 .* t .* t - cos(t);
               0 + t - t; % v
               cos(t);
               -t + sin(t);
               1 + t - t; % R
               0 + t - t;
               0 + t - t;
               0 + t - t;
               cos(t);
               sin(t);
               0 + t - t;
               -sin(t);
               cos(t);
               1 + t - t; % f
               0 + t - t; % p
               -sin(t);
               cos(t);
               0 + t - t; % psi_N
               0 + t - t; % psi_S
               0 + t - t; % dp
               -cos(t);
               -sin(t)];
% Initial state.
e3 = [0; 0; 1];
x0 = [0.5; 0; 0];
v0 = [0; 0; 0];
R0 = expm(pi / 3 * hat(e3));
f0 = 0.9;
X0 = [x0; v0; reshape(R0, 9, 1); f0];

p0 = f0 * R0 * e3;
if not(p0(1) == 0 && p0(2) == 0 && p0(3) <= 0)
    chart_pre = 'N';
else
    chart_pre = 'S';
end

tspan = [0, 7];
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

% Simulation.
[t, X] = ode45(@(t,X) sys(t,X,Xd_fun), tspan, X0, opts);

%% Plot.
close all;
Xd = Xd_fun(t')';
err = X - Xd(:, 1:16);

figure(1);
err_x = err(:, 1:3);
err_v = err(:, 4:6);
plot(t, sqrt(sum(err_x .* err_x, 2)));
hold on;
plot(t, sqrt(sum(err_v .* err_v, 2)), '--');

figure(2);
err_R = err(:, 7:15);
plot(t, sqrt(sum(err_R .* err_R, 2)));
