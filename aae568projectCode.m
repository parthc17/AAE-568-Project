% AAE 568 Project Script

% Part 1 - Obtain trajectory around asteroid

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
    xHat(:,i) = xBar + delX;

end

figure; plot3(y(:,1),y(:,2),y(:,3))
figure; plot3(xHat(1,:),xHat(2,:),xHat(3,:))
figure;
subplot(3,1,1); plot(t,y(:,1) - xHat(1,:)')
subplot(3,1,2); plot(t,y(:,2) - xHat(2,:)')
subplot(3,1,3); plot(t,y(:,3) - xHat(3,:)')

function dx = bennuProp(t,x,muBody)

    r = norm(x(1:3)); v = norm(x(4:6)); phi = reshape(x(7:42), [6 6]);
    radBennu = 250; % m
    j2 = 3.9257534110e-2;
    j3 = 1.4711698072e-2;

    % J2 Acceleration 
    j2leadingTerm = (-3*muBody*j2*radBennu^2) / (2*r^5);
    j2BCI = j2leadingTerm * [x(1)*(1-5*(x(3)^2/r^2)); x(2)*(1-5*(x(3)^2/r^2)); x(3)*(3-5*(x(3)^2/r^2))];

    % J3 Acceleration
    j3termXY = 5*muBody*radBennu^3*j3 / (2*r^7);
    j3termZ = muBody*radBennu^3*j3 / (2*r^5);
    j3BCI = [j3termXY * (7*x(3)^2 / r^2 - 3) * x(3)*x(1); j3termXY * (7*x(3)^2 / r^2 - 3) * x(3)*x(2); j3termZ * (3 - 30*x(3)^2/r^2 + 35*x(3)^4/r^4)];

    a = -muBody * x(1:3) / r^3 + j2BCI + j3BCI;

    A = zeros(6,6); A(1:3,4:6) = eye(3);
    dfdr = muBody * (3*x(1:3)*x(1:3)' / r^5 - eye(3)/r^3);
    dj2dx = -3*muBody*j2*radBennu^2 / (2*r^5) * [5*x(1)*(7*x(3)^2/r^2 - 1)/r^2 * x(1:3)' + (1-5*x(3)^2/r^2) * [1,0,0] - 10*x(1)*x(3)/r^2 * [0,0,1]];
    dj2dy = -3*muBody*j2*radBennu^2 / (2*r^5) * [5*x(2)*(7*x(3)^2/r^2 - 1)/r^2 * x(1:3)' + (1-5*x(3)^2/r^2) * [0,1,0] - 10*x(1)*x(3)/r^2 * [0,0,1]];
    dj2dz = -3*muBody*j2*radBennu^2 / (2*r^5) * [5*x(1)*(7*x(3)^2/r^2 - 3)/r^2 * x(1:3)' + 3*(1-5*x(3)^2/r^2) * [0,0,1]];
    A(4:6,1:3) = dfdr + [dj2dx;dj2dy;dj2dz];

    phiDot = A * phi;
    
    dx = [x(4:6); a; reshape(phiDot, [36 1])];

end