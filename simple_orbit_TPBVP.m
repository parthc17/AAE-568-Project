clearvars
close all

%% Setup/constants
muBennu = 5.2; % m^3/s^2
tol = 1e-12; opts = odeset('RelTol', tol, 'AbsTol', tol);
orbit = struct;
orbit.a = 1.5e3; orbit.e = 0.05; orbit.i = deg2rad(0); 
orbit.raan = 0; orbit.w = 0; orbit.f = 0;
[r0,v0] = keplerian2eci(orbit.a,orbit.e,orbit.i,orbit.raan,orbit.w,orbit.f,muBennu);
x0 = [r0; v0];
xf = [-r0; -v0];
period = 2*pi*sqrt(orbit.a^3/muBennu); % s
tf = period/2;
t = linspace(0, tf, 500);

l_char = 1000;
mu_char = muBennu;
t_char = sqrt(1/muBennu * l_char^3);
v_char = l_char/t_char;
a_char = l_char/t_char^2;

% NONDIM
x_scale = [l_char, l_char, l_char, v_char, v_char, v_char]';
x0_nd = x0 ./ x_scale;
xf_nd = xf ./ x_scale;
t_nd = t / t_char;
mu_nd = 1;

umax = 1;

%% Solve TPBVP
solinit = bvpinit(t_nd, @(t) guess_ytraj(t, x0_nd, period/t_char));  % initial guess p = 1
sol1 = bvp4c(@(t,x) state_coststate_dynamics(t,x, mu_nd, umax/a_char), ...
            @(y0, yf) bounds(y0, yf, x0_nd, xf_nd), solinit)
sol = bvp4c(@(t,x) state_coststate_dynamics(t,x, mu_nd, umax/a_char), ...
            @(y0, yf) bounds(y0, yf, x0_nd, xf_nd), sol1)
y_ref = sol.y(1:6,:) .* x_scale;
lm_ref = sol.y(7:12,:);
B = [zeros(3); eye(3)];
u_ref = min_energy_control(lm_ref, B, umax/a_char) * a_char;
y_ref(1:6,end) - x0;


[t,y2] = ode45(@(t,x) keplerian_dynamics(t, x, muBennu), t, [r0;v0], opts);

%% plot

%plot control
figure()
plot(t, u_ref)
xlabel("Time (s)"), ylabel("u [m/s2]"), grid(), legend(["ux", "uy", "uz"])
figure()
plot(t, lm_ref)
xlabel("Time (s)"), ylabel("Costate [ND]"), grid(),

% Plot 3d plot
y_g = guess_ytraj(t, x0, period);
figure; 
plot3(y_ref(1,:),y_ref(2,:),y_ref(3,:), DisplayName="TPBVP traj"); hold on;
plot3(x0(1), x0(2), x0(3), "or", DisplayName="x0")
plot3(xf(1), xf(2), xf(3), "*r", DisplayName="xf")
% plot3(y2(:,1),y2(:,2),y2(:,3), DisplayName="Uncontrolled More Perturbed");

fv = stlread('g_06290mm_spc_obj_0000n00000_v008.stl');
% Visualize correctly
pp = patch('Faces', fv.ConnectivityList, ...
          'Vertices', fv.Points, ...
          'FaceColor', [0.8 0.8 0.8], ...
          'EdgeColor', 'none', ...
          'FaceLighting', 'gouraud');

plot3(y_g(1,:),y_g(2,:),y_g(3,:), DisplayName="BVP Guess"); 
plot3(y2(:,1),y2(:,2),y2(:,3), ":", DisplayName="No control"); 
% plot3(y_ref(1,:),y_ref(2,:),y_ref(3,:), DisplayName="Controlled Reference"); 
legend();
hold off;

% Adjust view
axis equal
camlight headlight
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)'); title('Asteroid Orbit Control Truth'); 


%% Functions
function dXdt = keplerian_dynamics(t, X, mu) % not used?
    r = norm(X(1:3)); v = norm(X(4:6));
    a = -mu * X(1:3) / r^3; 
    dXdt = [X(4:6); a];
end
function u = min_energy_control(lm, B, umax) % u* optimal control input for min energy
    p = (-1/2*lm'*B)'; 
    u = p./vecnorm(p, 2, 1) .* clip(vecnorm(p, 2, 1), 0, umax);
    
end
function u = min_fuel_control(lm, B, umax) % u* optimal control input for min energy
    p = (-1/2*lm'*B)'; 
    S = (vecnorm(p, 2, 1) > 1);
    u = p./vecnorm(p, 2, 1) .* (S * umax);
end
function y = guess_ytraj(t, x0, period) % initial guess traj - rotate x0 to make a circle
    % global x0 tf
    % tf = t(end);
    n = length(t);
    w = 2*pi/period;
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
function dxdt = state_coststate_dynamics(t, x, mu, umax) % dynamics for TPBVP - mostly copied from Benuprop
    
    % if nondim
    %     lscale = 1/l_char;
    %     vscale = 1/(l_char/t_char);
    %     ascale = 1/(l_char/t_char^2);
    %     tscale = 1/t_char;
    % else
    %     lscale = 1;
    %     vscale = 1;
    %     ascale = 1;
    %     tscale = 1;
    % end

    % with J2
    % Constants
    r = norm(x(1:3)); v = norm(x(4:6));
    lm = x(7:12);

    ag = -mu * x(1:3) / r^3;%
    dfdr = mu * (3*x(1:3)*x(1:3)' / r^5 - eye(3)/r^3);

    % costate
    A = zeros(6,6); A(1:3,4:6) = eye(3);
    % A(4:6,1:3) = dfdr + dj2dr;
    A(4:6,1:3) = dfdr;
    dlmdt = - A' * lm;

    % control
    B = [zeros(3); eye(3)];

    u = min_energy_control(lm, B, umax);
    % u = min_fuel_control(lm, B, umax);

    dxdt = [x(4:6); ag + u; dlmdt];

end

function resid = bounds(y0, yf, x0, xf) % boundary conditions for TPBVP
    % global x0
    lmf = yf(7:12);
    % match xf to xo to make periodic orbit
    resid = [y0(1:6) - x0
             yf(1:6) - xf];
    a=0;
end