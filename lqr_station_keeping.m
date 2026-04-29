clearvars
clc
close all
%{
Aidan edits 4/16:
- Made sure it was being computed fully discretely
- Changed Ak and Bk so that they are just from tk to t(k+1) instead of t0 to tk
Hannah edits 4/23:
- cleaned up code
Hannah edits 4/27:
use min fuel
%}


%% Setup
rng(1)
sig_a_truth = 1e-6;   % Process Noise, m/s^2, tune this
muBennu = 5.2; % m^3/s^2
tol = 1e-12; opts = odeset('RelTol', tol, 'AbsTol', tol);
% orbit = struct;
% orbit.a = 1.5e3; orbit.e = 0.05; orbit.i = deg2rad(45); 
% orbit.raan = 0; orbit.w = 0; orbit.f = 0;
% [r0,v0] = keplerian2eci(orbit.a,orbit.e,orbit.i,orbit.raan,orbit.w,orbit.f,muBennu);
% x0 = [r0; v0];
% period = 2*pi*sqrt(orbit.a^3/muBennu); % s
% tf = 10*period;
% nt=5000;
% t = linspace(0, tf, nt);
nx = 6; nu = 3;
periods = 10;
tscale = 3600; 

%% Reference traj
% [t,y1] = ode45(@(t,x) bennuProp_old_ABC(t,x,muBennu), t, [r0;v0; reshape(eye(6),[nx^2 1]); reshape(zeros(nx, nu),[nx*nu 1])], opts); % target
ref = load("min_fuel_tpbvb.mat");
u_ref = ref.traj.u;
x_ref = ref.traj.x;
x0 = ref.traj.x0;
X = [x0; reshape(eye(6),[nx^2 1]); reshape(zeros(nx, nu),[nx*nu 1])];
y1 = zeros([60, size(ref.traj.x, 1)])';
t = ref.traj.t;
nt = size(t, 1);
np = nt;
for k = 1:(nt-1)
    % Do propagation - hold control constant over [t(k), t(k+1)]
    y1(k,:) = X;
    [~, Y] = ode45(@(tt,xx) bennuProp_old_ABC(tt, xx, muBennu, u_ref(:,k)), ...
                       [t(k) t(k+1)], X, opts);
    X = Y(end,:);
    X(1:6) = x_ref(k+1,:);
end
y1(nt,:) = X;
% x_ref = y1(:,1:6);
% y1 = repmat(y1(1:nt-1,:), periods, 1);
%% Get Ak, Bk, expand
phi_hist = y1(:,(nx+1):(nx+nx^2));
phi_hist = matrixify(phi_hist, 6, 6);
Bk_int = matrixify(y1(:,(nx+nx^2+1):(nx+nx^2+nx*nu)), nx, nu);
Bk_factor = pagemtimes(phi_hist, Bk_int);

% Make actual Ak and Bk for (tk, tk+1)
Ak_hist = zeros(size(phi_hist));
Bk_hist = zeros(size(Bk_factor));
% last item in hist will be zeros - discard
for k = 1:(size(phi_hist, 3)-1)
    Ak_hist(:,:,k) = phi_hist(:,:,k+1) / phi_hist(:,:,k);
    Bk_hist(:,:,k) = Bk_factor(:,:,k+1) - Ak_hist(:,:,k) * Bk_factor(:,:,k);
end

% repeat for multiple orbits
nt = (size(y1, 1)-1)*periods;
dtk = ref.traj.t(2)-ref.traj.t(1);
tf = periods * ref.traj.t(end);
t = linspace(0, tf, tf/dtk);

% repeat, but get rid of last one because it would duplicate IC and is 0
Ak_hist = repmat(Ak_hist(:,:,1:end-1), 1, 1, periods);
Bk_hist = repmat(Bk_hist(:,:,1:end-1), 1, 1, periods);
x_ref = repmat(x_ref(1:end-1, :), periods, 1);
u_ref = repmat(u_ref(:, 1:end-1), 1, periods);

% compute Ks for LQR (backwards in time)
% Q = diag([1,1,1,.5,.5,.5]) * 1;
% R = diag([1,1,1]) * 1;

max_pos_err = .5;      % m    - max acceptable position error
max_vel_err = 0.02;    % m/s  - max acceptable velocity error  
max_u       = 1e-8;    % m/s^2 - max acceptable control input

Q = diag(1 ./ [max_pos_err, max_pos_err, max_pos_err, ...
               max_vel_err, max_vel_err, max_vel_err].^2);
R = diag(1 ./ [max_u, max_u, max_u].^2);
Qf = Q; 

K_hist = zeros(nu, nx, nt);
Pk_hist = zeros(nx, nx, nt);
Pk = Qf;
Pk_hist(:,:,nt) = Qf;
for k = nt-1:-1:1
    % Ak = phi_hist(:,:,k+1) / phi_hist(:,:,k);
    % Bk = Bk_factor(:,:,k+1) - Ak * Bk_factor(:,:,k);
    Ak = Ak_hist(:,:,k);
    Bk = Bk_hist(:,:,k);

    S = R + Bk' * Pk * Bk;
    K_hist(:,:,k) = S \ (Bk' * Pk * Ak);

    Pk = Q + Ak' * Pk * Ak - Ak' * Pk * Bk * (S \ (Bk' * Pk * Ak));
    Pk_hist(:,:,k) = Pk;
end

% mpc stuff
N = floor(0.25*np)
% N = 5;

%% LQR Controlled 
% u_ref = zeros(nu, nt);

x_lqr = zeros(nx, nt);
x_mpc = zeros(nx, nt);
u_hist_lqr = zeros(nu, nt-1);
u_hist_mpc = zeros(nu, nt-1);

% Perturbed initial condition
% d0 = [0, 0, 0]';
% x_lqr(:,1) = [r0 + d0; v0];
% x_mpc(:,1) = [r0 + d0; v0];
x_lqr(:,1) = x0;
x_mpc(:,1) = x0;

% lqr
tic;
for k = 1:nt-1
    % tracking error at sample time tk
    % err_k = x_lqr(:,k) - x_ref(k,:)';

    % discrete-time TVLQR control
    u_k = u_ref(:,k) - K_hist(:,:,k) * (x_lqr(:,k) - x_ref(k,:)');
    u_hist_lqr(:,k) = u_k;
    
    % Do propagation - hold control constant over [t(k), t(k+1)]
    aNoise_k = sig_a_truth * randn(3,1);
    [~, y_seg_lqr] = ode45(@(tt,xx) bennuProp_controlled(tt, xx, muBennu, aNoise_k, u_k), ...
                       [t(k) t(k+1)], x_lqr(:,k), opts);

    x_lqr(:,k+1) = y_seg_lqr(end,:)';
end
fprintf("LQR Time: %.2f s\n", toc)

% mpc
tic;
for k = 1:nt-1
    % tracking error at sample time tk
    % err_k = x_lqr(:,k) - x_ref(k,:)';

    % MPC LQR
    P_N = Pk_hist(:,:,min(k+N, nt));
    K_mpc = mpc_control_linearized(Ak_hist, Bk_hist, k, R, Q, P_N, N, nt);
    u_k_mpc = u_ref(:,k) - K_mpc * (x_mpc(:,k) - x_ref(k,:)');
    u_hist_mpc(:,k) = u_k_mpc;
    
    % Do propagation - hold control constant over [t(k), t(k+1)]
    aNoise_k = sig_a_truth * randn(3,1);
    [~, y_seg_mpc] = ode45(@(tt,xx) bennuProp_controlled(tt, xx, muBennu, aNoise_k, u_k_mpc), ...
                   [t(k) t(k+1)], x_mpc(:,k), opts);

    x_mpc(:,k+1) = y_seg_mpc(end,:)';

end
fprintf("MPC Time: %.2f s\n", toc)


y3 = x_lqr.';
y4 = x_mpc.';

%% MPC (receding horizon LQR) controlled



%% Open Loop Traj
% w0 = randn([3,1]);
[t,y2] = ode45(@(t,x) bennuProp(t,x,muBennu,false), t, x0, opts); % open loop

err1 = y3(:,1:6) - x_ref;
err2 = y4(:,1:6) - x_ref;

%% Tracking error plots

figure('Name','Tracking Error vs Time');

subplot(2,1,1)
% plot(t, err1(:,1:3), 'LineWidth', 1.2)
plot(t/tscale, vecnorm(err1(:,1:3), 2, 2), 'LineWidth', 1.2), hold on
plot(t/tscale, vecnorm(err2(:,1:3), 2, 2), 'LineWidth', 1.2), hold off
grid on
xlabel('Time (hours)')
ylabel('Position Error (m)')
title('Position Tracking Error')
% legend('e_x', 'e_y', 'e_z', 'Location', 'best')
legend(["LQR", "MPC"])

subplot(2,1,2)
% plot(t, err1(:,4:6), 'LineWidth', 1.2)
plot(t/tscale, vecnorm(err1(:,4:6), 2, 2), 'LineWidth', 1.2), hold on
plot(t/tscale, vecnorm(err2(:,4:6), 2, 2), 'LineWidth', 1.2), hold off
grid on
xlabel('Time (s)')
ylabel('Velocity Error (m/s)')
title('Velocity Tracking Error')
% legend('e_{v_x}', 'e_{v_y}', 'e_{v_z}', 'Location', 'best')
legend(["LQR", "MPC"])
fig = gcf();
saveas(fig, "Controller_error.png")

%% Control history plot
% u_hist_lqr = u_ref - squeeze(pagemtimes(K_hist, reshape(err1', nx, 1, nt)));

figure('Name','LQR Control Correction vs Time');
plot(t(1:end-1)/tscale, squeeze(vecnorm(u_hist_lqr, 2, 1)), 'LineWidth', 0.8), hold on
plot(t(1:end-1)/tscale, squeeze(vecnorm(u_hist_mpc, 2, 1)), 'LineWidth', 0.8), hold off
grid on
xlabel('Time (hours)')
ylabel('Control Acceleration (m/s^2)')
title('LQR Control Correction History')
% legend('u_x', 'u_y', 'u_z', 'Location', 'best')
legend(["LQR", "MPC"])

fig = gcf();
saveas(fig, "Controller_u.png")

%% 3D trajectory comparison
figure('Name','3D Trajectory Comparison');
plot3(y1(:,1), y1(:,2), y1(:,3), 'LineWidth', 1.5, ...
    'DisplayName', 'Reference: Uncontrolled J2 Perturbed');
hold on
plot3(y2(:,1), y2(:,2), y2(:,3), 'LineWidth', 1.5, ...
    'DisplayName', 'Uncontrolled Perturbed');
plot3(y3(:,1), y3(:,2), y3(:,3), '--', 'LineWidth', 1.5, ...
    'DisplayName', 'LQR Controlled');
plot3(y4(:,1), y4(:,2), y4(:,3), '--', 'LineWidth', 1.5, ...
    'DisplayName', 'MPC Controlled');

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

err_norm = vecnorm((y3(:,1:6) - x_ref)',2,1);
figure; plot(t, err_norm); grid on
title('Tracking Error Norm');
xlabel('Time (s)');
ylabel('||x_{ctrl} - x_{ref}||');


%% functions
function dx = bennuProp_old_ABC(t,x,muBody, u)

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
    % j3termXY = 5*muBody*radBennu^3*j3 / (2*r^7);
    % j3termZ = muBody*radBennu^3*j3 / (2*r^5);
    % j3BCI = [j3termXY * (7*x(3)^2 / r^2 - 3) * x(3)*x(1); j3termXY * (7*x(3)^2 / r^2 - 3) * x(3)*x(2); j3termZ * (3 - 30*x(3)^2/r^2 + 35*x(3)^4/r^4)];

    a = -muBody * x(1:3) / r^3 + j2BCI + u; % + j3BCI;

    A = zeros(6,6); A(1:3,4:6) = eye(3);
    A(4:6,1:3) = dfdr + dj2dr;

    B = [zeros(3); eye(3)];

    phiDot = A * phi; % Ak
    Bkdot = phi \ B;
    
    dx = [x(4:6); a; reshape(phiDot, [36 1]); reshape(Bkdot, [18, 1])];

end

function K = mpc_control_linearized(Ak_hist, Bk_hist, k, R, Q, P_N, N, kmax)
    % Make big matrices
    n = size(Ak_hist, 1);
    m = size(Bk_hist, 2);
    Qb = kron(eye(N), Q); % N blkdiag of Q
    Qb(end-n+1:end, end-n+1:end) = Qb(end-n+1:end, end-n+1:end) + P_N;
    Rb = kron(eye(N), R); % N blkdiag of R
    Kn = [eye(m), zeros(m, m*(N-1))];
    H = zeros(N*n, n);
    G = zeros(n*N, m*N);
    % G_row = [B, zeros(n, m*(N-1))];
    % G_row = 
    for i = 1:min(N,kmax-k) % rows
        A_stack = eye(n);
        G_row = zeros(n, m*(N));
        for j = i:-1:1 % nonzero G cols or H rows
            % B = Bk_hist(k+i-1);
            G_row(:, j:(j+m-1)) = A_stack * Bk_hist(:,:,k+j-1);
            if j ~= 1 % the last update of A_stack isn't used so don't compute it so it doesnt break the end
                A_stack = A_stack * Ak_hist(:,:,k+j);
            end
            % H(((i-1)*n+1):(i*n), :) = A^i;
        end
        H(((i-1)*n+1):(i*n), :) = A_stack * Ak_hist(:,:,k);
        G(((i-1)*n+1):(i*n), :) = G_row;
            % G_row = [A*G_row(:, 1:m), G_row(:,1:m*(N-1))];
    end
    
    F = G' * Qb * H;
    L = G' * Qb * G + Rb;
    K = Kn * (L \ F);
end


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
