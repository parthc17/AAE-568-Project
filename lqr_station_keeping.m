clearvars
clc
close all
%{
Aidan edits 4/16:
- Made sure it was being computed fully discretely
- Changed Ak and Bk so that they are just from tk to t(k+1) instead of t0 to tk
%}


%% Setup
rng(1)
sig_a_truth = 1e-7;   % Process Noise, m/s^2, tune this
muBennu = 5.2; % m^3/s^2
tol = 1e-12; opts = odeset('RelTol', tol, 'AbsTol', tol);
orbit = struct;
orbit.a = 1.5e3; orbit.e = 0.05; orbit.i = deg2rad(45); 
orbit.raan = 0; orbit.w = 0; orbit.f = 0;
[r0,v0] = keplerian2eci(orbit.a,orbit.e,orbit.i,orbit.raan,orbit.w,orbit.f,muBennu);
x0 = [r0; v0];
period = 2*pi*sqrt(orbit.a^3/muBennu); % s
tf = 10*period;
nt=5000;
t = linspace(0, tf, nt);
dtk = tf/nt;
nx = 6; nu = 3;

%% Reference traj
[t,y1] = ode45(@(t,x) bennuProp_old_ABC(t,x,muBennu), t, [r0;v0; reshape(eye(6),[nx^2 1]); reshape(zeros(nx, nu),[nx*nu 1])], opts); % target
x_ref = y1(:,1:6);
phi_hist = y1(:,(nx+1):(nx+nx^2));
phi_hist = matrixify(phi_hist, 6, 6);
Bk_int = matrixify(y1(:,(nx+nx^2+1):(nx+nx^2+nx*nu)), nx, nu);
Bk_hist = pagemtimes(phi_hist, Bk_int);

% compute Ks (backwards in time)
% k_out = ode45(@(t,x) k_ode(t,x, A, B, P, Q, R), fliplr(tspan), K0, options);
Q = diag([1,1,1,.5,.5,.5]) * 1;
R = diag([1,1,1]) * 1;

% K_hist = zeros(nu, nx, nt);
% Pk = Q; % final Pk
% for k=fliplr(1:nt) % solve backwards
%     Ak = phi_hist(:,:,k);
%     Bk = Bk_hist(:,:,k);
%     K_hist(:,:,k) = (R + Bk'*Pk*Bk) \ Bk'*Pk*Ak;
%     Pk = Q + Ak'*Pk*Ak - (Ak'*Pk*Bk) / (R + Bk'*Pk*Bk) * (Bk'*Pk*Ak);
% end

Qf = Q*2; 
K_hist = zeros(nu, nx, nt);
Pk = Qf;
for k = nt-1:-1:1
    Ak = phi_hist(:,:,k+1) / phi_hist(:,:,k);
    Bk = Bk_hist(:,:,k+1) - Ak * Bk_hist(:,:,k);

    S = R + Bk' * Pk * Bk;
    K_hist(:,:,k) = S \ (Bk' * Pk * Ak);

    Pk = Q + Ak' * Pk * Ak - Ak' * Pk * Bk * (S \ (Bk' * Pk * Ak));
end

%% Controlled 
u_ref = zeros(nu, nt);

x_ctrl = zeros(nx, nt);
u_hist = zeros(nu, nt-1);

% Perturbed initial condition
x_ctrl(:,1) = [r0 + [10;10;10]; v0];

for k = 1:nt-1
    % tracking error at sample time tk
    err_k = x_ctrl(:,k) - x_ref(k,:)';

    % discrete-time TVLQR control
    u_k = u_ref(:,k) - K_hist(:,:,k) * err_k;
    u_hist(:,k) = u_k;
    
    % hold control constant over [t(k), t(k+1)]
    aNoise_k = sig_a_truth * randn(3,1);
    [~, y_seg] = ode45(@(tt,xx) bennuProp_controlled(tt, xx, muBennu, aNoise_k, u_k), ...
                       [t(k) t(k+1)], x_ctrl(:,k), opts);

    x_ctrl(:,k+1) = y_seg(end,:)';
end

y3 = x_ctrl.';

%% Open Loop Traj
[t,y2] = ode45(@(t,x) bennuProp(t,x,muBennu,false), t, [r0 + [10;10;10];v0], opts); % open loop

err = y3(:,1:6) - y1(:,1:6);

%% Tracking error plots
err = y3(:,1:6) - y1(:,1:6);

figure('Name','Tracking Error vs Time');

subplot(2,1,1)
plot(t, err(:,1:3), 'LineWidth', 1.2)
grid on
xlabel('Time (s)')
ylabel('Position Error (m)')
title('Position Tracking Error')
legend('e_x', 'e_y', 'e_z', 'Location', 'best')

subplot(2,1,2)
plot(t, err(:,4:6), 'LineWidth', 1.2)
grid on
xlabel('Time (s)')
ylabel('Velocity Error (m/s)')
title('Velocity Tracking Error')
legend('e_{v_x}', 'e_{v_y}', 'e_{v_z}', 'Location', 'best')


%% Control history plot
u_hist = u_ref - squeeze(pagemtimes(K_hist, reshape(err', nx, 1, nt)));

figure('Name','LQR Control Correction vs Time');
plot(t, squeeze(u_hist), 'LineWidth', 1.2)
grid on
xlabel('Time (s)')
ylabel('Control Acceleration (m/s^2)')
title('LQR Control Correction History')
legend('u_x', 'u_y', 'u_z', 'Location', 'best')


%% 3D trajectory comparison
figure('Name','3D Trajectory Comparison');
plot3(y1(:,1), y1(:,2), y1(:,3), 'LineWidth', 1.5, ...
    'DisplayName', 'Reference: Uncontrolled J2 Perturbed');
hold on
plot3(y2(:,1), y2(:,2), y2(:,3), 'LineWidth', 1.5, ...
    'DisplayName', 'Uncontrolled Perturbed');
plot3(y3(:,1), y3(:,2), y3(:,3), '--', 'LineWidth', 1.5, ...
    'DisplayName', 'LQR Controlled');

fv = stlread('g_06290mm_spc_obj_0000n00000_v008.stl');
patch('Faces', fv.ConnectivityList, ...
      'Vertices', fv.Points, ...
      'FaceColor', [0.8 0.8 0.8], ...
      'EdgeColor', 'none', ...
      'FaceLighting', 'gouraud', ...
      'DisplayName', 'Bennu Shape');

grid on
axis equal
camlight headlight
xlabel('X Position (m)')
ylabel('Y Position (m)')
zlabel('Z Position (m)')
title('Asteroid Orbit Control Trajectories')
legend('Location', 'best')
hold off

err_norm = vecnorm((y3(:,1:6) - y1(:,1:6))',2,1);
figure; plot(t, err_norm); grid on
title('Tracking Error Norm');
xlabel('Time (s)');
ylabel('||x_{ctrl} - x_{ref}||');


%% dynamics
function dx = bennuProp_old_ABC(t,x,muBody)

    % Constants
    r = norm(x(1:3)); v = norm(x(4:6)); phi = reshape(x(7:42), [6 6]);
    radBennu = 250; % m
    j2 = 3.9257534110e-2;
    j3 = 1.4711698072e-2;

    % J2 Acceleration and Partials
    j2leadingTerm = (-3*muBody*j2*radBennu^2) / (2*r^5);
    j2BCI = j2leadingTerm * [x(1)*(1-5*(x(3)^2/r^2)); x(2)*(1-5*(x(3)^2/r^2)); x(3)*(3-5*(x(3)^2/r^2))];

    dfdr = muBody * (3*x(1:3)*x(1:3)' / r^5 - eye(3)/r^3);
    dj2dx = -3*muBody*j2*radBennu^2 / (2*r^5) * [5*x(1)*(7*x(3)^2/r^2 - 1)/r^2 * x(1:3)' + (1-5*x(3)^2/r^2) * [1,0,0] - 10*x(1)*x(3)/r^2 * [0,0,1]];
    dj2dy = -3*muBody*j2*radBennu^2 / (2*r^5) * [5*x(2)*(7*x(3)^2/r^2 - 1)/r^2 * x(1:3)' + (1-5*x(3)^2/r^2) * [0,1,0] - 10*x(1)*x(3)/r^2 * [0,0,1]];
    dj2dz = -3*muBody*j2*radBennu^2 / (2*r^5) * [5*x(1)*(7*x(3)^2/r^2 - 3)/r^2 * x(1:3)' + 3*(1-5*x(3)^2/r^2) * [0,0,1]];
    dj2dr = [dj2dx;dj2dy;dj2dz];

    % J3 Acceleration and Partials
    j3termXY = 5*muBody*radBennu^3*j3 / (2*r^7);
    j3termZ = muBody*radBennu^3*j3 / (2*r^5);
    j3BCI = [j3termXY * (7*x(3)^2 / r^2 - 3) * x(3)*x(1); j3termXY * (7*x(3)^2 / r^2 - 3) * x(3)*x(2); j3termZ * (3 - 30*x(3)^2/r^2 + 35*x(3)^4/r^4)];

    a = -muBody * x(1:3) / r^3;% + j2BCI;% + j3BCI;

    A = zeros(6,6); A(1:3,4:6) = eye(3);
    A(4:6,1:3) = dfdr;% + dj2dr;

    B = [zeros(3); eye(3)];

    phiDot = A * phi; % Ak
    Bkdot = phi \ B;
    
    dx = [x(4:6); a; reshape(phiDot, [36 1]); reshape(Bkdot, [18, 1])];

end

% function dxdt_lqr = bennuProp_lqr(t, x, muBody, x_ref, u_ref, K_hist, dtk)
%     % x_ref, u_ref shape [n x nk] (transpose from ode)
%     % NOTE: Using K from Ak, Bk from REFERENCE, not actual traj. 
%     % but this could be modified to calculate Ak, Bk on actual traj
%     nt = size(x_ref, 2); % number of k steps
%     tspan = linspace(0,nt*dtk, nt);
%     k = min(floor(t/dtk) + 1, nt);
% 
%     dxdt = bennuProp(t,x, muBody, false); % get dynamics from propagator function
%     err = x - interp1(tspan, x_ref', t, "cubic")';
%     % NOTE: should K be interpolated? Or maybe just use the continuous one
%     % instead
%     uk = u_ref(:,k) - K_hist(:,:,k) * err;
%     B = [zeros(3); eye(3)]; % B =/= Bk!
%     dxdt_lqr = dxdt + B * uk;
% end

%% Stuff from ACA
% x0 = [r0;v0;log(m0)];
% Y0 = [x0;phi0];
% traj = ode45(@(t,x) eomstm(t,x, g, alpha), tspan, Y0, options);
% Y = deval(traj, 10);
% Ak = matrixify(Y(nx+1:nx+nx^2), nx, nx)
% Bk = Ak * matrixify(Y(nx+nx^2+1:nx+nx^2+nx*nu), nx, nu)
% Ck = Ak * matrixify(Y(nx+nx^2+nx*nu+1:nx+nx^2+nx*nu+nx), nx, 1)

% Numerical Solution
% k_out = ode45(@(t,x) k_ode(t,x, A, B, P, Q, R), fliplr(tspan), K0, options);
% 
% K_hist = flip(matrixify(k_out.y', n, n),3); % flipped so in time order
% 
% tk = k_out.x;
% [t1,x] = ode45(@(t,x) x_ode(t,x, A, B, P, Q, R, k_out), flip(tk), x0, options);
% x1 = reshape(x',n,1,[]);
% lm1 = pagemtimes(K_hist, x1);
% % u = -R\(P'*x + B'*y);
% u1 = pagemtimes(-inv(R), (pagemtimes(P, 'transpose', x1, 'none') + pagemtimes(B, 'transpose', lm1, 'none')));

% function k_dot = k_ode(t,x, A, B, P, Q, R)
%     % x = [x, K]'
%     n = length(A);
%     % m = length(R);
%     % x = X(1:n);
%     K = reshape(x, n, n);
% 
%     % u = -R\(P'+B'*K)*x; % omptimal control u
%     % x_dot = A*x + B*u;
% 
%     k_dot = -K * (A - B/R*(P' + B'*K)) - Q + P/R*(P' + B'*K) - A'*K;
%     k_dot = reshape(k_dot, 16, 1);
%     % dXdt = [x_dot; k_dot];
% end
% 
% function x_dot = x_ode(t,x, A, B, P, Q, R, K_hist)
%     % x = [x, K]'
%     n = length(A);
%     % m = length(R);
%     % x = X(1:n);
%     K = deval(K_hist, t);
%     K = reshape(K, n, n);
% 
%     u = -R\(P'+B'*K)*x; % omptimal control u
%     x_dot = A*x + B*u;
% 
%     % k_dot = -K * (A - B/R*(P' + B'*K)) - Q + P/R*(P' + B'*K) - A'*K;
%     % k_dot = reshape(k_dot, 16, 1);
%     % dXdt = [x_dot; k_dot];
% end

%Dynamics
function dx = bennuProp(t,x,muBody,propSTM,aNoise)
    if nargin < 5
        aNoise = [0;0;0];
    elseif nargin < 4
        propSTM = false;
    end

    % Constants
    r = norm(x(1:3)); v = norm(x(4:6)); 
    radBennu = 250; % m
    j2 = 3.9257534110e-2;
    % j3 = 1.4711698072e-2;

    % J2 Acceleration and Partials
    j2leadingTerm = (-3*muBody*j2*radBennu^2) / (2*r^5);
    j2BCI = j2leadingTerm * [x(1)*(1-5*(x(3)^2/r^2)); x(2)*(1-5*(x(3)^2/r^2)); x(3)*(3-5*(x(3)^2/r^2))];

    dfdr = muBody * (3*x(1:3)*x(1:3)' / r^5 - eye(3)/r^3);
    dj2dx = -3*muBody*j2*radBennu^2 / (2*r^5) * [5*x(1)*(7*x(3)^2/r^2 - 1)/r^2 * x(1:3)' + (1-5*x(3)^2/r^2) * [1,0,0] - 10*x(1)*x(3)/r^2 * [0,0,1]];
    dj2dy = -3*muBody*j2*radBennu^2 / (2*r^5) * [5*x(2)*(7*x(3)^2/r^2 - 1)/r^2 * x(1:3)' + (1-5*x(3)^2/r^2) * [0,1,0] - 10*x(1)*x(3)/r^2 * [0,0,1]];
    dj2dz = -3*muBody*j2*radBennu^2 / (2*r^5) * [5*x(1)*(7*x(3)^2/r^2 - 3)/r^2 * x(1:3)' + 3*(1-5*x(3)^2/r^2) * [0,0,1]];
    dj2dr = [dj2dx;dj2dy;dj2dz];

    % % J3 Acceleration and Partials
    % j3termXY = 5*muBody*radBennu^3*j3 / (2*r^7);
    % j3termZ = muBody*radBennu^3*j3 / (2*r^5);
    % j3BCI = [j3termXY * (7*x(3)^2 / r^2 - 3) * x(3)*x(1); j3termXY * (7*x(3)^2 / r^2 - 3) * x(3)*x(2); j3termZ * (3 - 30*x(3)^2/r^2 + 35*x(3)^4/r^4)];

    % Perturbation due to Sun's gravity
    G = 6.674e-11; % Universal gravitational constant, N*m^2/kg^2
    m2 = 1.989e30; % Mass of the Sun, kg
    mInAU = 1.496192602435979E+11; % m/AU
    r_ben2sun = [1.36 * mInAU; 0; 0];  % Bennu-to-Sun vector, m
    r_sc = x(1:3);  %Spacecraft position relative to Bennu, m
    d = r_ben2sun - r_sc;   % Sun-relative vector from spacecraft to Sun, m
    a_sun = G*m2 * ( d / norm(d)^3 - r_ben2sun / norm(r_ben2sun)^3 ); % Third-body perturbation acceleration due to Sun

    % Partial derivative of Sun third-body perturbation wrt spacecraft position
    daSun_dr = G*m2 * ( 3*(d*d')/norm(d)^5 - eye(3)/norm(d)^3 );

    %Perturbation due to Solar Pressure
    SF = 1353; %Solar radiation constant, W/m^2
    c = 3*10^8; %Speed of light, m/s
    P_SR_1AU = SF/c; %Solar Pressure per m^2, at 1 AU
    P_SR = P_SR_1AU * (mInAU/norm(r_ben2sun))^2; %Adjustment for being further than 1 AU
    c_R = 0.6; %Reflectivity of the satellite
    A_exposed = 14; %Exposed surface area of the satellite to the sun, m^2
    m = 800; %Mass of the satellite, kg
    a_SR = -P_SR*c_R*A_exposed*d/(m*norm(d)); %Perturbation due to solar radiation pressure

    %Partial derivative of SRP perturbation
    daSRP_dr = -P_SR*c_R*A_exposed/m*(d*d'/norm(d)^3 - eye(3)/norm(d));

    % a = -muBody * x(1:3) / r^3 + j2BCI + a_sun + a_SR + aNoise; %All perturbations
    a = -muBody * x(1:3) / r^3 + j2BCI+ a_sun + a_SR + aNoise; %All perturbations

    if propSTM
        phi = reshape(x(7:42), [6 6]);
        A = zeros(6,6); A(1:3,4:6) = eye(3);
        A(4:6,1:3) = dfdr+ dj2dr + daSun_dr + daSRP_dr;% + daSRP_dr;  %State matrix with perturbations

        phiDot = A * phi; %Derivative of STM
    
        dx = [x(4:6); a; reshape(phiDot, [36 1])];
    else
        dx = [x(4:6); a];
    end

end

function dx = bennuProp_controlled(t, x, muBody, aNoise_k, u)
    dx_nom = bennuProp(t, x, muBody, false, aNoise_k);  % nominal nonlinear dynamics
    B = [zeros(3,3); eye(3)];
    dx = dx_nom + B*u;
end
