clearvars
close all

%% Setup
muBennu = 5.2; % m^3/s^2
tol = 1e-12; opts = odeset('RelTol', tol, 'AbsTol', tol);
orbit = struct;
orbit.a = 1.5e3; orbit.e = 0.05; orbit.i = deg2rad(45); 
orbit.raan = 0; orbit.w = 0; orbit.f = 0;
[r0,v0] = keplerian2eci(orbit.a,orbit.e,orbit.i,orbit.raan,orbit.w,orbit.f,muBennu);
x0 = [r0; v0]
period = 2*pi*sqrt(orbit.a^3/muBennu); % s
tf = period;
nt=500;
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
Q = eye(nx);
R = eye(nu);

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
    Ak = phi_hist(:,:,k);
    Bk = Bk_hist(:,:,k);

    S = R + Bk' * Pk * Bk;
    K_hist(:,:,k) = S \ (Bk' * Pk * Ak);

    Pk = Q + Ak' * Pk * Ak - Ak' * Pk * Bk * (S \ (Bk' * Pk * Ak));
end

%% Controlled 
u_ref = zeros(nu, nt); % should come from OL control BV4PC
[t, y3] = ode45(@(t,x) bennuProp_lqr(t, x, muBennu, x_ref', u_ref, K_hist, dtk), ...
        t, [r0;v0], opts); % open loop

%% Open Loop Traj
[t,y2] = ode45(@(t,x) bennuProp(t,x,muBennu,false), t, [r0;v0], opts); % open loop

err = y3(:,1:6) - y1(:,1:6);
figure();
subplot(2,1,1)
plot(t, err(:,1:3)), grid(), ylabel("Position Error [m]")
subplot(2,1,2)
plot(t, err(:,4:6)), grid(), ylabel("Velocity Error [m/s]")

u_hist = u_ref - squeeze(pagemtimes(K_hist(), reshape(err', nx, 1, nt)));
figure();
plot(t, squeeze(u_hist))
grid(), ylabel("LQR Control Correction [m/s2]")


%% Plot results
% Plot True Dynamics
figure; 
plot3(y1(:,1),y1(:,2),y1(:,3), DisplayName="Uncontrolled J2 Perturbed (Ref)"); hold on;
plot3(y2(:,1),y2(:,2),y2(:,3), DisplayName="Uncontrolled More Perturbed");
plot3(y3(:,1),y3(:,2),y3(:,3), "--", DisplayName="LQR Controlled");

fv = stlread('g_06290mm_spc_obj_0000n00000_v008.stl');
% Visualize correctly
pp = patch('Faces', fv.ConnectivityList, ...
          'Vertices', fv.Points, ...
          'FaceColor', [0.8 0.8 0.8], ...
          'EdgeColor', 'none', ...
          'FaceLighting', 'gouraud');

% plot3(y_g(1,:),y_g(2,:),y_g(3,:), DisplayName="BVP Guess"); 
% plot3(y_ref(1,:),y_ref(2,:),y_ref(3,:), DisplayName="Controlled Reference"); 
legend();
hold off;

% Adjust view
axis equal
camlight headlight
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)'); title('Asteroid Orbit Control Truth'); 

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

    a = -muBody * x(1:3) / r^3 + j2BCI;% + j3BCI;

    A = zeros(6,6); A(1:3,4:6) = eye(3);
    A(4:6,1:3) = dfdr + dj2dr;

    B = [zeros(3); eye(3)];

    phiDot = A * phi; % Ak
    Bkdot = phi \ B;
    
    dx = [x(4:6); a; reshape(phiDot, [36 1]); reshape(Bkdot, [18, 1])];

end

function dxdt_lqr = bennuProp_lqr(t, x, muBody, x_ref, u_ref, K_hist, dtk)
    % x_ref, u_ref shape [n x nk] (transpose from ode)
    % NOTE: Using K from Ak, Bk from REFERENCE, not actual traj. 
    % but this could be modified to calculate Ak, Bk on actual traj
    nt = size(x_ref, 2); % number of k steps
    tspan = linspace(0,nt*dtk, nt);
    k = min(floor(t/dtk) + 1, nt);

    dxdt = bennuProp(t,x, muBody, false); % get dynamics from propagator function
    err = x - interp1(tspan, x_ref', t, "cubic")';
    % NOTE: should K be interpolated? Or maybe just use the continuous one
    % instead
    uk = u_ref(:,k) - K_hist(:,:,k) * err;
    B = [zeros(3); eye(3)]; % B =/= Bk!
    dxdt_lqr = dxdt + B * uk;
end

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

