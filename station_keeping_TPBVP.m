clearvars
close all;
%% Setup
global x0 tf umax t_char l_char v_char a_char nondim
nondim = true;
% umax = 1
umax = 1e-6;
muBennu = 5.2; % m^3/s^2
    radBennu = 250; % m

tol = 1e-12; opts = odeset('RelTol', tol, 'AbsTol', tol);
orbit = struct;
orbit.a = 1.5e3; orbit.e = 0.05; orbit.i = deg2rad(45); 
orbit.raan = 0; orbit.w = 0; orbit.f = 0;
[r0,v0] = keplerian2eci(orbit.a,orbit.e,orbit.i,orbit.raan,orbit.w,orbit.f,muBennu);
x0 = [r0; v0];
period = 2*pi*sqrt(orbit.a^3/muBennu); % s
tf = period;
t = linspace(0, tf, 500);
dt = tf/500;
B = [zeros(3); eye(3)];

l_char = 1000;
mu_char = muBennu;
mu_nd = 1;
t_char = sqrt(1/muBennu * l_char^3);
v_char = l_char/t_char;
a_char = l_char/t_char^2;

% NONDIM
x_scale = [l_char, l_char, l_char, v_char, v_char, v_char]';
x0_nd = x0 ./ x_scale;
t_nd = t / t_char;

%% Solve TPBVP
solinit = bvpinit(t_nd, @(t) guess_ytraj(t, x0_nd, t_nd(end)));  % initial guess p = 1
% do min energy
sol = bvp4c(@(t,x) state_coststate_dynamics_energy(t,x, umax/a_char, mu_nd, radBennu/l_char), ...
            @(y0, yf) bounds(y0, yf, x0_nd), solinit)
sol = bvp4c(@(t,x) state_coststate_dynamics_energy(t,x, umax/a_char, mu_nd, radBennu/l_char), ...
            @(y0, yf) bounds(y0, yf, x0_nd), sol)
lm_ref = sol.y(7:12,:);
u_ref = min_energy_control(lm_ref, B, umax/a_char) * a_char;
dv = sum(vecnorm(u_ref,2,1))*dt;
fprintf("Min energy dv = %.6f m/s\n", dv)
% do min fuel
rhos = [0.5, 0.1, 0.01, 0.001];
for p=rhos
    fprintf("Solveing rho = %.2f\n", p)
    sol = bvp4c(@(t,x) state_coststate_dynamics_fuel(t,x, umax/a_char, mu_nd, radBennu/l_char, p), ...
                @(y0, yf) bounds(y0, yf, x0_nd), sol)
end
sol = bvp4c(@(t,x) state_coststate_dynamics_fuel(t,x, umax/a_char, mu_nd, radBennu/l_char, 0.1), ...
            @(y0, yf) bounds(y0, yf, x0_nd), sol)
%
y_ref = sol.y(1:6,:) .* x_scale;
lm_ref = sol.y(7:12,:);
u_ref = min_fuel_control(lm_ref, B, umax/a_char, p) * a_char;
% u_ref = min_energy_control(lm_ref, B, umax/a_char) * a_char;
dv = sum(vecnorm(u_ref,2,1))*dt;
fprintf("Min fuel dv = %.6f m/s\n", dv)

y_ref(1:6,end) - x0


% Truth Propagation
[t,y1] = ode45(@(t,x) bennuProp_old(t,x,muBennu), t, [r0;v0; reshape(eye(6),[36 1])], opts);
y1(end,1:6)' - x0;

[t,y2] = ode45(@(t,x) bennuProp(t,x,muBennu,true), t, [r0;v0; reshape(eye(6),[36 1])], opts);

% guess prop
y_g = guess_ytraj(t, x0, tf);

%% plot control
figure()
plot(t, u_ref), hold on
plot(t, vecnorm(u_ref, 2, 1), "k")
plot(xlim(), [umax, umax], "r--"), hold off
xlabel("Time (s)"), ylabel("u [m/s2]"), grid(), legend(["ux", "uy", "uz"])
figure()
plot(t, lm_ref)
xlabel("Time (s)"), ylabel("Costate [ND]"), grid(),

%% Plot 3d plot
figure; 
plot3(y1(:,1),y1(:,2),y1(:,3), DisplayName="Uncontrolled J2 Perturbed"); hold on;
plot3(y2(:,1),y2(:,2),y2(:,3), DisplayName="Uncontrolled More Perturbed");

fv = stlread('g_06290mm_spc_obj_0000n00000_v008.stl');
% Visualize correctly
pp = patch('Faces', fv.ConnectivityList, ...
          'Vertices', fv.Points, ...
          'FaceColor', [0.8 0.8 0.8], ...
          'EdgeColor', 'none', ...
          'FaceLighting', 'gouraud');

plot3(y_g(1,:),y_g(2,:),y_g(3,:), DisplayName="BVP Guess"); 
plot3(y_ref(1,:),y_ref(2,:),y_ref(3,:), DisplayName="Controlled Reference"); 
legend();
hold off;

% Adjust view
axis equal
camlight headlight
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)'); title('Asteroid Orbit Control Truth'); 

%% Functions
function dx = bennuProp_old(t,x,muBody) % with less perturbations

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
function dXdt = keplerian_dynamics(t, X) % not used?
    r = norm(X(1:3)); v = norm(X(4:6));
    a = -muBody * X(1:3) / r^3; 
    dXdt = [X(4:6); a];
end

function [a, A] = bennu_dynamics(x, muBody, radBennu)
    r = norm(x(1:3)); v = norm(x(4:6));
    j2 = 3.9257534110e-2;

    ag = -muBody * x(1:3) / r^3;%
    dfdr = muBody * (3*x(1:3)*x(1:3)' / r^5 - eye(3)/r^3);

    % J2 Acceleration and Partials
    j2leadingTerm = (-3*muBody*j2*radBennu^2) / (2*r^5);
    j2BCI = j2leadingTerm * [x(1)*(1-5*(x(3)^2/r^2)); x(2)*(1-5*(x(3)^2/r^2)); x(3)*(3-5*(x(3)^2/r^2))];

    dj2dx = -3*muBody*j2*radBennu^2 / (2*r^5) * [5*x(1)*(7*x(3)^2/r^2 - 1)/r^2 * x(1:3)' + (1-5*x(3)^2/r^2) * [1,0,0] - 10*x(1)*x(3)/r^2 * [0,0,1]];
    dj2dy = -3*muBody*j2*radBennu^2 / (2*r^5) * [5*x(2)*(7*x(3)^2/r^2 - 1)/r^2 * x(1:3)' + (1-5*x(3)^2/r^2) * [0,1,0] - 10*x(1)*x(3)/r^2 * [0,0,1]];
    dj2dz = -3*muBody*j2*radBennu^2 / (2*r^5) * [5*x(1)*(7*x(3)^2/r^2 - 3)/r^2 * x(1:3)' + 3*(1-5*x(3)^2/r^2) * [0,0,1]];
    dj2dr = [dj2dx;dj2dy;dj2dz];

    A = zeros(6,6); A(1:3,4:6) = eye(3);
    A(4:6,1:3) = dfdr + dj2dr;
    a = ag + j2BCI;
end

function dxdt = state_coststate_dynamics_energy(t, x, umax, muBody, radBennu) % dynamics for TPBVP - mostly copied from Benuprop    
    % with J2
    lm = x(7:12);

    [a, A] = bennu_dynamics(x, muBody, radBennu);
    
    % costate
    dlmdt = - A' * lm;

    % control
    B = [zeros(3); eye(3)];

    u = min_energy_control(lm, B, umax);
    % u = min_fuel_control(lm, B, umax);

    dxdt = [x(4:6); a + u; dlmdt];
end

function dxdt = state_coststate_dynamics_fuel(t, x, umax, muBody, radBennu, rho) % dynamics for TPBVP - mostly copied from Benuprop    
    % with J2
    lm = x(7:12);

    [a, A] = bennu_dynamics(x, muBody, radBennu);
    
    % costate
    dlmdt = - A' * lm;

    % control
    B = [zeros(3); eye(3)];

    % u = min_energy_control(lm, B, umax);
    u = min_fuel_control(lm, B, umax, rho);

    dxdt = [x(4:6); a + u; dlmdt];
end

function resid = bounds(y0, yf, x0) % boundary conditions for TPBVP
    % global x0
    lmf = yf(7:12);
    % match xf to xo to make periodic orbit
    resid = [y0(1:6) - x0
             yf(1:6) - x0];
    a=0;
end

function y = guess_ytraj(t, x0, tf) % initial guess traj - rotate x0 to make a circle
    % global x0 tf
    % tf = t(end);
    n = length(t);
    w = 2*pi/tf;
    angles = -w*t;
    r = norm(x0(1:3));
    ax = cross(x0(4:6), x0(1:3));
    y = ones(12, n)*-.1;
    for i=1:n
        y(1:3,i) = rodrigues_rot(x0(1:3), ax, angles(i));
        y(4:6,i) = rodrigues_rot(x0(4:6), ax, angles(i));
    end
    % y = [x0 + (xf - x0)*t
    %     1
    %     w*t];
end

function u = min_energy_control(lm, B, umax) % u* optimal control input for min energy
    p = (-1/2*lm'*B)'; 
    u = p./vecnorm(p, 2, 1) .* clip(vecnorm(p, 2, 1), 0, umax);
    
end

function u = min_fuel_control(lm, B, umax, rho) % u* optimal control input for min energy
    p = (-1/2*lm'*B)'; 
    % S = (vecnorm(p, 2, 1) > 1);
    S = vecnorm(p, 2, 1)-1;
    u = p./vecnorm(p, 2, 1) .* (umax/2 * (1 + tanh(S/rho)));
end


