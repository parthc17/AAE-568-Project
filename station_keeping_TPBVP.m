function dx = bennuProp_old(t,x,muBody)

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
function dXdt = keplerian_dynamics(t, X)
    r = norm(X(1:3)); v = norm(X(4:6));
    a = -muBody * X(1:3) / r^3; 
    dXdt = [X(4:6); a];
end

function dxdt = state_coststate_dynamics(t, x, umax)
    global l_char t_char nondim
    % nondim = True;
    
    if nondim
        lscale = 1/l_char;
        vscale = 1/(l_char/t_char);
        ascale = 1/(l_char/t_char^2);
        tscale = 1/t_char;
    else
        lscale = 1;
        vscale = 1;
        ascale = 1;
        tscale = 1;
    end
    % with J2
    % Constants
    muBody = 5.2 * (lscale^3/tscale^2); % m^3/s^2
    r = norm(x(1:3)); v = norm(x(4:6));
    lm = x(7:12);
    radBennu = 250 * lscale; % m
    j2 = 3.9257534110e-2; % ND

    ag = -muBody * x(1:3) / r^3;%
    dfdr = muBody * (3*x(1:3)*x(1:3)' / r^5 - eye(3)/r^3);

    % J2 Acceleration and Partials
    j2leadingTerm = (-3*muBody*j2*radBennu^2) / (2*r^5);
    j2BCI = j2leadingTerm * [x(1)*(1-5*(x(3)^2/r^2)); x(2)*(1-5*(x(3)^2/r^2)); x(3)*(3-5*(x(3)^2/r^2))];

    dj2dx = -3*muBody*j2*radBennu^2 / (2*r^5) * [5*x(1)*(7*x(3)^2/r^2 - 1)/r^2 * x(1:3)' + (1-5*x(3)^2/r^2) * [1,0,0] - 10*x(1)*x(3)/r^2 * [0,0,1]];
    dj2dy = -3*muBody*j2*radBennu^2 / (2*r^5) * [5*x(2)*(7*x(3)^2/r^2 - 1)/r^2 * x(1:3)' + (1-5*x(3)^2/r^2) * [0,1,0] - 10*x(1)*x(3)/r^2 * [0,0,1]];
    dj2dz = -3*muBody*j2*radBennu^2 / (2*r^5) * [5*x(1)*(7*x(3)^2/r^2 - 3)/r^2 * x(1:3)' + 3*(1-5*x(3)^2/r^2) * [0,0,1]];
    dj2dr = [dj2dx;dj2dy;dj2dz];

    % costate
    A = zeros(6,6); A(1:3,4:6) = eye(3);
    A(4:6,1:3) = dfdr + dj2dr;
    dlmdt = - A' * lm;

    % control
    B = [zeros(3); eye(3)];

    u = min_energy_control(lm, B, umax);
    % u = min_fuel_control(lm, B, umax);

    dxdt = [x(4:6); ag + j2BCI + u; dlmdt];
end

function resid = bounds(y0, yf, x0)
    % global x0
    lmf = yf(7:12);
    resid = [y0(1:6) - x0
             yf(1:6) - x0];
    a=0;
end

function y = guess_ytraj(t, x0, tf)
    % global x0 tf
    % tf = t(end);
    n = length(t);
    w = 2*pi/tf;
    angles = w*t;
    r = norm(x0(1:3));
    ax = cross(x0(4:6), x0(1:3));
    y = ones(12, n)*.1;
    for i=1:n
        y(1:3,i) = rodrigues_rot(x0(1:3), ax, angles(i));
        y(4:6,i) = rodrigues_rot(x0(4:6), ax, angles(i));
    end
    % y = [x0 + (xf - x0)*t
    %     1
    %     w*t];
end

function u = min_energy_control(lm, B, umax)
    p = (-1/2*lm'*B)'; 
    u = p./vecnorm(p, 2, 1) .* clip(vecnorm(p, 2, 1), 0, umax);
    
end

function u = min_fuel_control(lm, B, umax)
    p = (-1/2*lm'*B)'; 
    S = (vecnorm(p, 2, 1) > 1);
    u = p./vecnorm(p, 2, 1) .* (S * umax);
end

global x0 tf umax t_char l_char v_char a_char nondim
nondim = true;
umax = 1e-5
% umax = 0.01
muBennu = 5.2; % m^3/s^2
tol = 1e-12; opts = odeset('RelTol', tol, 'AbsTol', tol);
orbit = struct;
orbit.a = 1.5e3; orbit.e = 0.05; orbit.i = deg2rad(45); 
orbit.raan = 0; orbit.w = 0; orbit.f = 0;
[r0,v0] = keplerian2eci(orbit.a,orbit.e,orbit.i,orbit.raan,orbit.w,orbit.f,muBennu);
x0 = [r0; v0]
period = 2*pi*sqrt(orbit.a^3/muBennu); % s
tf = period
t = linspace(0, tf, 500);

l_char = 1000
t_char = 3600
v_char = l_char/t_char;
a_char = l_char/t_char^2;

% NONDIM
x_scale = [l_char, l_char, l_char, v_char, v_char, v_char]';
x0_nd = x0 ./ x_scale;
t_nd = t / t_char;

% TPBVP
solinit = bvpinit(t_nd, @(t) guess_ytraj(t, x0_nd, t_nd(end)));  % initial guess p = 1
sol = bvp4c(@(t,x) state_coststate_dynamics(t,x, umax/a_char), ...
            @(y0, yf) bounds(y0, yf, x0_nd), solinit)
y_ref = sol.y(1:6,:) .* x_scale;
lm_ref = sol.y(7:12,:)
% p_ref = -1/2*y_ref(10:12,:);
% u_ref = p_ref/vecnorm(p_ref, 2, 1) * clip(vecnorm(p_ref, 2, 1), -umax, umax);
B = [zeros(3); eye(3)];
u_ref = min_energy_control(lm_ref, B, umax) * a_char;
y_ref(1:6,end) - x0

% Truth Propagation
[t,y1] = ode45(@(t,x) bennuProp_old(t,x,muBennu), t, [r0;v0; reshape(eye(6),[36 1])], opts);
y1(end,1:6)' - x0

[t,y2] = ode45(@(t,x) bennuProp(t,x,muBennu), t, [r0;v0; reshape(eye(6),[36 1])], opts);

% guess prop
y_g = guess_ytraj(t, x0, tf);

% plot control
figure()
plot(t, u_ref)
xlabel("Time (s)"), ylabel("u [m/s2]"), grid(), legend(["ux", "uy", "uz"])
figure()
plot(t, lm_ref)
xlabel("Time (s)"), ylabel("Costate [ND]"), grid(),

% Plot True Dynamics
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

