function dX = sys(t, X, Xd)
global chart_pre
% Constants.
m = 1;
g = 1;
e3 = [0; 0; 1];

% State.
x = X(1:3);
v = X(4:6);
R = reshape(X(7:15), 3, 3);
f = X(16);

% Desired state.
Xd_t = Xd(t);
xd = Xd_t(1:3);
vd = Xd_t(4:6);
% Rd = reshape(Xd_t(7:15), 3, 3);
% fd = Xd_t(16);

% Compute eta.
pd = Xd_t(17:19);
psi_N_d = Xd_t(20);
psi_S_d = Xd_t(21);
dpd = Xd_t(22:24);

A = f * R;
p = A * e3;
p_norm = norm(p);
K = blkdiag(2, eye(3));
w = dpd - 30 * (x - xd) - 28 * (v - vd) - 9 * (p - pd);

if chart_pre == 'N'
    if sqrt(p(1)*p(1)+p(2)*p(2)) < p_norm * sin(pi/6) && p(3) < 0
        eta = US(A, p, w, psi_S_d, K);
        chart_pre = 'S';
        disp(['Switch from UN to US at time ', num2str(t)]);
    else
        eta = UN(A, p, w, psi_N_d, K);
    end
else
    if sqrt(p(1)*p(1)+p(2)*p(2)) < p_norm * sin(pi/6) && p(3) > 0
        eta = UN(A, p, w, psi_N_d, K);
        chart_pre = 'N';
        disp(['Switch from US to UN at time ', num2str(t)]);
    else
        eta = US(A, p, w, psi_S_d, K);
    end
end

% Control.
u0 = eta(1);
omega = eta(2:4);

% Dynamics.
dx = v;
dv = -g * e3 + f / m * R * e3;
dR = hat(omega) * R;
df = f * u0;

dX = [dx; dv; reshape(dR,9,1); df];
end

function eta = US(A, p, w, psi_d, K)
psi_S = atan2(-A(1,2)-A(2,1), A(1,1)-A(2,2));
w0_S = -3 * (psi_S - psi_d);

p_norm = norm(p);
FS = [   0, p(1)/(2*p_norm-p(3)), p(2)/(2*(p_norm-p(3))),  -0.5;
      p(1),                    0,                   p(3), -p(2);
      p(2),                -p(3),                      0,  p(1);
      p(3),                 p(2),                  -p(1),     0];

eta = pinv(K * FS) * [w0_S; w];
end

function eta = UN(A, p, w, psi_d, K)
psi_N = atan2(A(2,1)-A(1,2), A(1,1)+A(2,2));
w0_N = -3 * (psi_N - psi_d);

p_norm = norm(p);
FN = [   0, p(1)/(2*p_norm+p(3)), p(2)/(2*(p_norm+p(3))),   0.5;
      p(1),                    0,                   p(3), -p(2);
      p(2),                -p(3),                      0,  p(1);
      p(3),                 p(2),                  -p(1),     0];

eta = pinv(K * FN) * [w0_N; w];
end