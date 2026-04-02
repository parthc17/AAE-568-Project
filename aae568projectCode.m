% AAE 568 Project Script

% Part 1 - Obtain trajectory to Bennu - work in progress

% earthState0 = struct;
% earthState0.a = 1; % au
% earthState0.e = 0.0169;
% earthState0.i = 4.001e-5;
% earthState0.raan = 2.4473;
% earthState0.w = 5.6647;
% earthState0.f = 1.8858;
% 
% bennuState0 = struct;
% bennuState0.a = 1.078; % au
% bennuState0.e = 0.8270;
% bennuState0.i = 0.398;
% bennuState0.raan = 1.5356;
% bennuState0.w = 0.5489;
% bennuState0.f = 3.8412;
% 
% muSun = 1.32712440018e11;
% kmInAU = 1.496192602435979E+08; % km/AU
% 
% % Starred values
% lStar = kmInAU; % lStar = 1.496e+8; % km In AU
% tStar = sqrt(lStar^3/muSunDim); % s
% 
% muSun = muSun / (lStar^3/tStar^2);
% 
% % Time span
% t0 = 0; 
% tf = 8; % Minimum Fuel Time
% timesp = linspace(t0, tf*tStar, 10000);
% tol = 1e-12; opts = odeset('RelTol', tol, 'AbsTol', tol);
% 
% [rBennu0,vBennu0] = keplerian2eci(bennuState0.a,bennuState0.e,bennuState0.i,bennuState0.raan,bennuState0.w,bennuState0.f,muSun);
% [tBennu, xBennuFinal] = ode45(@(t,x) cartesian(t,x,muSun), timesp, [rBennu0;vBennu0], opts); % Obtain final Bennu state
% 



%% Part 2 - Orbit Determination for Bennu orbiting

% Constants
muBennu = 5.2; % m^3/s^2
tol = 1e-12; opts = odeset('RelTol', tol, 'AbsTol', tol);
ixx = 1.8130e9; 
iyy = 1.8836e9;
izz = 2.0334e9; % kg km^2
mass = 7.7e10;

j2 = 3.9257534110e-2;
j3 = 1.4711698072e-2;
j4 = 3.0760445246e-2;

% Initial Orbit
orbit = struct;
orbit.a = 1.5e3; orbit.e = 0.05; orbit.i = deg2rad(45); 
orbit.raan = 0; orbit.w = 0; orbit.f = 0;

[r0,v0] = keplerian2eci(orbit.a,orbit.e,orbit.i,orbit.raan,orbit.w,orbit.f,muBennu);
cov0 = diag([10,10,10,0.1,0.1,0.1]);

% Time
time = linspace(0,50*86400,5000);
period = 2*pi*sqrt(orbit.a^3/muBennu) / 86400

% Truth Propagation
[t,y] = ode45(@(t,x) bennuProp(t,x,muBennu), time, [r0;v0; reshape(eye(6),[36 1])], opts);

% Noisy Propagation
[tN,yN] = ode45(@(t,x) bennuProp(t,x,muBennu), time, [r0+randn(3,1);v0; reshape(eye(6),[36 1])], opts);

% EKF

for i = 1:length(t)
    phi(:,:,i) = reshape(y(i,7:42),[6 6]);
end

H = [eye(3) zeros(3,3)];
Q = eye(6) / 10; Q(4:6,4:6) = Q(4:6,4:6) / 10;
R = 0.1*eye(3);
Parr = zeros(6,6,length(t)); Parr(:,:,1) = cov0;
for i = 1:length(t)
    P = Parr(:,:,i);
    phi = reshape(y(i,7:42),[6 6]);
    xBar = yN(i,1:6)';
    covBar = phi*P*phi' + Q;

    meas = y(i,1:3)' + (randn(3,1) / 10);
    % meas(1:3) = meas(1:3) + (randn(3,1) / 10);
    % meas(4:6) = meas(4:6) + (randn(3,1) / 1000);

    b = meas - H*xBar;
    K = covBar * H' * inv(H*covBar*H' + R);
    delX = K*b;
    Phat(:,:,i) = covBar - K*H*covBar;
    Pxx(i) = Phat(1,1,i);
    Pyy(i) = Phat(2,2,i);
    Pzz(i) = Phat(3,3,i);
    xHat(:,i) = xBar + delX;

end

% Plot True Dynamics
figure; plot3(y(:,1),y(:,2),y(:,3)); hold on;
fv = stlread('g_06290mm_spc_obj_0000n00000_v008.stl');
% Visualize correctly
pp = patch('Faces', fv.ConnectivityList, ...
          'Vertices', fv.Points, ...
          'FaceColor', [0.8 0.8 0.8], ...
          'EdgeColor', 'none', ...
          'FaceLighting', 'gouraud');
% Adjust view
axis equal
camlight headlight
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)'); title('Asteroid Orbit Control'); 

% Plot EKF Results
figure; plot3(xHat(1,:),xHat(2,:),xHat(3,:)); hold on;
fv = stlread('g_06290mm_spc_obj_0000n00000_v008.stl');
% Visualize correctly
pp = patch('Faces', fv.ConnectivityList, ...
          'Vertices', fv.Points, ...
          'FaceColor', [0.8 0.8 0.8], ...
          'EdgeColor', 'none', ...
          'FaceLighting', 'gouraud');
% Adjust view
axis equal
camlight headlight
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)'); title('Asteroid Orbit Control'); 

% Plot Differences between Truth and EKF
figure;
subplot(3,1,1); plot(t,y(:,1) - xHat(1,:)')
subplot(3,1,2); plot(t,y(:,2) - xHat(2,:)')
subplot(3,1,3); plot(t,y(:,3) - xHat(3,:)')

% Plot sigma bounds
figure;
subplot(3,1,1); plot(t,xHat(1,:)); hold on; plot(t, xHat(1,:) + 3 * sqrt(Pxx)); plot(t, xHat(1,:) - 3 * sqrt(Pxx))
subplot(3,1,2); plot(t,xHat(2,:)); hold on; plot(t, xHat(2,:) + 3 * sqrt(Pyy)); plot(t, xHat(2,:) - 3 * sqrt(Pyy))
subplot(3,1,3); plot(t,xHat(3,:)); hold on; plot(t, xHat(3,:) + 3 * sqrt(Pzz)); plot(t, xHat(3,:) - 3 * sqrt(Pzz))

%% Transfer



%% Functions

function dx = bennuProp(t,x,muBody)

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

    phiDot = A * phi;
    
    dx = [x(4:6); a; reshape(phiDot, [36 1])];

end

function motion = BVP_ode(t, x, rho, uMax, imp, g, A, G0, muSun, thirdBodyStruct)

    % Set state and costate
    state = x(1:7);
    lambda = x(8:end);

    % Set state and costate vectors
    r = state(1:3); v = state(4:6); m = state(7);
    lr = lambda(1:3); lv = lambda(4:6); lm = lambda(7);

    if m < 0
        return;
    end

    % Control input setup
    uHatStar = -lv / norm(lv,2);
    S = 1 + lv' * uHatStar / m - lm / (imp*g);
    gammaStar = 0.5 * uMax * (1 + tanh(-S / rho));
    u = gammaStar * uHatStar;

    % dadx
    dadr  = -muSun * (eye(3) / norm(r)^3 - 3 * (r * r') / norm(r)^5);
    dadm = -u / m^2;

    % Acceleration
    accel = -r/(norm(r)^3) + u/m;

    % Derivative Values for integration
    xDot = [v; accel; -norm(u,2) / (imp*g)];
    dfdx = [zeros(3,3), eye(3), zeros(3,1);
            dadr, zeros(3,3), dadm;
            zeros(1,3), zeros(1,3), 0];
    lamDot = (-lambda'*dfdx)';

    motion = [xDot; lamDot];
end

function psi = BVP_BC(ya,yb)

    psi = [yb(1:6) - target;
           yb(14)];

end
