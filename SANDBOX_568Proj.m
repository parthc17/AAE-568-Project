clc; clear; close all;
% AAE 568 Project Script
%{
AIDAN EDITS 4/13:
- Removed the overwriting of the BC for Earth and Bennu locations within
the BC function
- Added function for initial guess trajectory instead of initial guess
state
- Decreased mesh size to help convergence
- Added a small floor to uHatStar
- Added small floor to mass
- Added function for using earth trajectory as initial guess
- Going to try removing mass constraint to get it to converge
- Removed terminal velocity constraints to help get it to converge
- Increase max thrust to help it converge
- Added terminal velocity back in
- Added mass back in
%}



% Part 1 - Obtain trajectory to Bennu - work in progress

muSun = 1.32712440018e11;
kmInAU = 1.496192602435979E+08; % km/AU

% Starred values
lStar = kmInAU; % lStar = 1.496e+8; % km In AU
tStar = sqrt(lStar^3/muSun); % s

earthState0 = struct;
earthState0.a = 1.496657326987069E+08 / lStar; % au
earthState0.e = 1.704313732350883E-02; earthState0.i = deg2rad(6.198205899446798E-03); earthState0.raan = deg2rad(1.799051298362160E+02);
earthState0.w = deg2rad(2.816178320319530E+02); earthState0.f = deg2rad(2.669012326892521E+02);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bennuState0 = struct;
bennuState0.a = 1.684375473676060E+08 / lStar; % au
bennuState0.e = 2.037604127989827E-01; bennuState0.i = deg2rad(6.032225719717332E+00); bennuState0.raan = deg2rad(1.949096214163919E+00);
bennuState0.w = deg2rad(6.640855196853940E+01); bennuState0.f = deg2rad(2.693898971751065E+02);

g0 = (9.80665 / 1000) / (lStar / tStar^2); % m/s^2 => nondim
isp = 4000 / (tStar); % s => nondim | 4000 s works
muSun = muSun / (lStar^3/tStar^2); % nondim mu Sun

% Time span
t0 = 0; 
tf = 8; % Minimum Fuel Time
%% CHANGE #1: timsesp changed to 80 increments. Basically making the mesh coarse
timesp = linspace(t0, tf, 80); % nondim time
tol = 1e-12; opts = odeset('RelTol', tol, 'AbsTol', tol);

[rBennu0,vBennu0] = keplerian2eci(bennuState0.a,bennuState0.e,bennuState0.i,bennuState0.raan,bennuState0.w,bennuState0.f,muSun);
% ^ is different from fsolve bc i used a different f value for some reason
[tBennu, xBennuFinal] = ode45(@(t,x) cartesian(t,x,muSun), timesp, [rBennu0;vBennu0], opts); % Obtain final Bennu state

%% CHANGE #2: bennuStateFinal is the end state of the bvp4c
bennuStateFinal = xBennuFinal(end,:);

%% CHANGE 3: xEarth is the initial state of the bvp4c
[rEarth0,vEarth0] = keplerian2eci(earthState0.a,earthState0.e,earthState0.i,earthState0.raan,earthState0.w,earthState0.f,muSun);
[tEarth, xEarth] = ode45(@(t,x) cartesian(t,x,muSun), timesp, [rEarth0;vEarth0], opts); % Obtain Earth Traj



options = bvpset('Stats','on','RelTol',1e-1);

%% CHANGE 4: This function generates a straight line trajectory from initial to final state
%guessFun = @(t) initialGuess(t, t0, tf, rEarth0, vEarth0, bennuStateFinal');

%% CHANGE 5: This function generates trajectory that is just the earth's orbit
 guessFun = @(t) earthTrajectoryGuess(t, tEarth, xEarth);

 %If you use a guess function, you still need to initialize it like a guess
 %state. This is the step where it takes that function and builds the
 %trajectory from it.
 %solinit = bvpinit(timesp, guessFun);
 %sol = bvp4c(@(t, x) BVP_ode(t, x, rho), @(ya, yb) BVP_BC(ya,yb,rEarth0,vEarth0,bennuStateFinal), solinit, options);

%% Once you get a converged solution, you can use it instead of a guess function
%% No need to initialize the BVP this way, the sol is already in correct format
% Initial_Guess = load('568ProjEarthToBennuGuess2_MaxInput061_TermVelIncluded.mat');
% sol = Initial_Guess.sol; 
%% By increasing rho and working our way down, we make it easier to converge. Large rho makes our control smoother, which helps with convergence
%% Once we get it working for smooth control, we use that solution for our next iteration, where we allow control to be more aggressive, and so on
% for rho = 20:-1:1
%     sol =  bvp4c(@(t, x) BVP_ode(t, x, rho), @(ya, yb) BVP_BC(ya,yb,rEarth0,vEarth0,bennuStateFinal), sol, options);
% end
% y = sol.y; 
% plot3(y(1,:),y(2,:),y(3,:));


%% Here is where I add mass back in. I am using the constant mass solution as my initial guess
% This means that I need take the constant mass solution, and add in a
% guess for how the mass should change over the trajectory
Initial_Guess = load('568ProjEarthToBennuGuess3_MaxInput061_MassNOTIncluded.mat');
solConstMass = Initial_Guess.sol;

m0 = 20; %Initial mass
mf_guess = 19.5;   %Guess for final mass
uMax_mass = 0.58;  %I moved uMax outside the function so it can be an input

%This takes the constant mass solution and adds the extra states for mass
%with a guess for them
massGuessFun = @(t) massGuess(t, solConstMass, m0, mf_guess);
solinit_mass = bvpinit(solConstMass.x, massGuessFun); 

% Start with easier/smoother rho first
rhoVals = [50 20 10 5 2 1];

%Get a initial converged solution before making control less smooth
sol_mass = bvp4c(@(t,x) BVP_ode_mass(t,x,rhoVals(1),uMax_mass), ...
                 @(ya,yb) BVP_BC_mass(ya,yb,rEarth0,vEarth0,bennuStateFinal,m0), ...
                 solinit_mass, options);

%Inch rho down to 1. Honestly might not be necessary anymore, I was more
%helpful for bringing residuals down significantly before I got a converged
%solution
for k = 2:length(rhoVals)
    rho = rhoVals(k);
    sol_mass = bvp4c(@(t,x) BVP_ode_mass(t,x,rho,uMax_mass), ...
                     @(ya,yb) BVP_BC_mass(ya,yb,rEarth0,vEarth0,bennuStateFinal,m0), ...
                     sol_mass, options);
end

%% Plot result
figure;
plot3(sol_mass.y(1,:), sol_mass.y(2,:), sol_mass.y(3,:), 'LineWidth', 1.5); hold on
plot3(rEarth0(1),rEarth0(2),rEarth0(3),'o','MarkerSize',7)
plot3(bennuStateFinal(1),bennuStateFinal(2),bennuStateFinal(3),'rx','MarkerSize',7)
plot3(xBennuFinal(:,1),xBennuFinal(:,2),xBennuFinal(:,3))
plot3(xEarth(:,1),xEarth(:,2),xEarth(:,3))
plot3(0,0,0,'go','MarkerSize',15)
axis equal
grid on
legend('Mass-inclusive transfer','Initial Position','Final Bennu Position','Bennu Orbit','Earth Orbit','Sun')
title('Earth-to-Bennu Transfer with Mass Dynamics')

%% Plot mass history just to check that its working
t_mass = sol_mass.x;
y_mass = sol_mass.y;

m_hist = y_mass(7,:); 

figure;
plot(t_mass, m_hist, 'LineWidth', 1.5);
grid on
xlabel('Time (nondimensional)');
ylabel('Mass (nondimensional)');
title('Spacecraft Mass vs Time');

%% Plot mass consumed
m0 = m_hist(1);
mass_used = m0 - m_hist;

figure;
plot(t_mass, mass_used, 'LineWidth', 1.5);
grid on
xlabel('Time (nondimensional)');
ylabel('Mass Consumed (nondimensional)');
title('Cumulative Mass Consumption vs Time');

%% Print quick diagnostics
fprintf('Initial mass: %.6f\n', m_hist(1));
fprintf('Final mass:   %.6f\n', m_hist(end));
fprintf('Mass used:    %.6f\n', m_hist(1) - m_hist(end));
fprintf('Minimum mass over trajectory: %.6f\n', min(m_hist));



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

%% CONSTANT MASS FUNCTIONS
function motion = BVP_ode_CONSTANTMASS(t, x, rho)

    g0 = 1.654184883569765e+03; % nondim
    isp = 7.962226975300293e-04; % nondim
    uMax = 0.61;
    muSun = 1;

    % Constant mass for simplified problem
    m = 20;

    % State and costate
    state = x(1:6);
    lambda = x(7:12);

    r = state(1:3);
    v = state(4:6);

    lr = lambda(1:3);
    lv = lambda(4:6);

    % Control input setup
    lvNorm = max(norm(lv,2), 1e-8);
    uHatStar = -lv / lvNorm;

    % Simplified switching function with constant mass
    S = 1 + lv' * uHatStar / m;
    gammaStar = 0.5 * uMax * (1 + tanh(-S / rho));
    u = gammaStar * uHatStar;

    % State dynamics
    dadr  = -muSun * (eye(3) / norm(r)^3 - 3 * (r * r') / norm(r)^5);
    accel = -r/(norm(r)^3) + u/m;

    xDot = [v;
            accel];

    % Jacobian of state dynamics wrt state [r;v]
    dfdx = [zeros(3,3), eye(3);
            dadr,       zeros(3,3)];

    % Costate dynamics
    lamDot = (-lambda' * dfdx)';

    motion = [xDot; lamDot];
end

function psi = BVP_BC_CONSTANTMASS(ya,yb,rEarth0,vEarth0, bennuStateFinal)
    psi = [ya(1:3) - rEarth0;
       ya(4:6) - vEarth0;
       yb(1:6) - bennuStateFinal'];

end

%% TRAJECTORY GUESSING FUNCTIONS
function z = initialGuess(t, t0, tf, r0, v0, xF)
    %Draws a straight line from state A to state B
    %Not necessarily physically realistic, but better than guessing a
    %single point as the solution trajectory

    s = (t - t0)/(tf - t0); %Basically tau from our hw, its time nondim...
    ...time for 0 to 1

    rF = xF(1:3); %Final position
    vF = xF(4:6); %Final vel

    rGuess = (1-s)*r0 + s*rF; %Straight line from position A to position B
    vGuess = (1-s)*v0 + s*vF; %Straight line from vel A to vel B

    lamGuess = zeros(6,1);
    lamGuess(4:6) = [-0.01; -0.01; 0];   %idk just small numbers tbh, arbitrary

    z = [rGuess; vGuess; lamGuess];
end

function z = earthTrajectoryGuess(t, tEarth, xEarth)
    % Interpolate Earth ballistic trajectory at time t
    % This states our initial trajectory as an orbit shape, which should
    % help us converge
    rGuess = interp1(tEarth, xEarth(:,1:3), t, 'linear')';
    vGuess = interp1(tEarth, xEarth(:,4:6), t, 'linear')';

    % Simple mass guess: slowly decreasing from initial mass
    s = (t - tEarth(1)) / (tEarth(end) - tEarth(1));

    % Small costate guess
    lamGuess = zeros(6,1);
    lamGuess(4:6) = [-0.01; -0.01; 0]; %Again, arbitrary, just want them small

    z = [rGuess; vGuess; lamGuess];
end

%% VARIABLE MASS FUNCTIONS
function motion = BVP_ode_mass(t, x, rho, uMax)

    g0 = 1.654184883569765e+03;   % nondim
    isp = 7.962226975300293e-04;  % nondim
    muSun = 1;

    % State and costate
    state = x(1:7);
    lambda = x(8:14);

    r = state(1:3);
    v = state(4:6);
    m = state(7);

    lr = lambda(1:3);
    lv = lambda(4:6);
    lm = lambda(7);

    % Numerical floors
    lvNorm = max(norm(lv,2), 1e-8);
    mEff   = max(m, 1e-6);

    % Control
    uHatStar  = -lv / lvNorm;
    S         = 1 + lv' * uHatStar / mEff - lm / (isp*g0);
    gammaStar = 0.5 * uMax * (1 + tanh(-S / rho));
    u         = gammaStar * uHatStar;

    % Dynamics
    dadr  = -muSun * (eye(3) / norm(r)^3 - 3 * (r * r') / norm(r)^5);
    dadm  = -u / mEff^2;
    accel = -r/(norm(r)^3) + u/mEff;

    xDot = [v;
            accel;
            -norm(u,2)/(isp*g0)];

    dfdx = [zeros(3,3), eye(3), zeros(3,1);
            dadr,       zeros(3,3), dadm;
            zeros(1,3), zeros(1,3), 0];

    lamDot = (-lambda' * dfdx)';

    motion = [xDot; lamDot];
end

function psi = BVP_BC_mass(ya, yb, rEarth0, vEarth0, bennuStateFinal, m0)

    psi = [ya(1:3) - rEarth0;
           ya(4:6) - vEarth0;
           ya(7)   - m0;
           yb(1:6) - bennuStateFinal';
           yb(14)];
end

function z = massGuess(t, solConstMass, m0, mf_guess)

    % Evaluate the converged constant-mass solution at time t
    y12 = deval(solConstMass, t);   % 12x1 vector

    % Normalized time
    s = (t - solConstMass.x(1)) / (solConstMass.x(end) - solConstMass.x(1));

    % Guessed mass history
    mGuess = m0 + (mf_guess - m0) * s;

    % Guessed lambda_m
    lmGuess = 0;

    % Build 14x1 vector:
    % [r(3); v(3); m; lambda_r(3); lambda_v(3); lambda_m]
    z = [y12(1:6);
         mGuess;
         y12(7:12);
         lmGuess];
end