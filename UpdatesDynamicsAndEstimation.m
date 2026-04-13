clc; clear; close all;
rng(1) %Set rng seed for repeatability
% AAE 568 Project Script
%{
---------------------------------------------------------------------------
AIDAN EDITS 4/10-4/11: 
- Added rng(1) at the top to make sure the randomness was the same every run
- Added titles/axis titles/legends to plots
- Combined the two trajectory plots into 1 so it could be seen how close
the two actually are
- Changed time to be constant time-step
- Added perturbations for sun's gravity and solar radiation pressure
- Modified EKF algorithm to more closely allign with algorithm from paper. 
Estimation is based on measurements and predicted trajectory, which is now
used to calculate the next predicted state. Previously the predicted state
was propogated independently of estimations and measurements.
- Modified error tuning matrices Q and R to place greater confidence in
predicted state over measured state.
- Added more noise to measurements
- Added process noise to truth trajectory to simulate unmodeled/unknown
perturbations in acceleration
- Satellite only takes measurements once every 5 minutes
AIDAN EDITS 4/12: 
- EKF Moved into a function so that other estimators can be added
- Added a secondary estimator that takes multiple measurements each time,
and then uses weighted batch least squares estimator to improve the
measurement, which is then used in the EKF
- Increased noise 
- Changed up some plots
-------------------------------------------------------------------------
%}

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

% Time (Constant time-step size)
dt = 100; %sec
time = 0:dt:(21.55*86400); %Roughly time until impact, uncontrolled
period = 2*pi*sqrt(orbit.a^3/muBennu) / 86400;

%% Truth Propagation------------------------------------------------------
sig_a_truth = 1e-7;   % Process Noise, m/s^2, tune this
yTruth = zeros(length(time),42);
yTruth(1,:) = [r0; v0; reshape(eye(6),[36 1])]';

for k = 1:length(time)-1
    aNoise_k = sig_a_truth * randn(3,1);

    [~, ySeg] = ode45(@(tt,xx) bennuProp(tt,xx,muBennu,aNoise_k), ...
                      [time(k) time(k+1)], yTruth(k,:)', opts);
    %Note: Each time step sees a random unknown perturbation that the
    %estimator needs to cope with. By adding this, our estimator's dynamics
    %is simulated as being imperfect, which makes the truth-based
    %measurements more important.
    yTruth(k+1,:) = ySeg(end,:);
end

t = time(:);
y = yTruth;

%% Estimation Setup------------------------------------------------------
sig = 2; %Standard Dev of Random noise, arbitrarily selected

H = [eye(3) zeros(3,3)]; %We can measure x, y, z
Q = diag([1e-6, 1e-6, 1e-6, 1e-10, 1e-10, 1e-10]); %Tune This, Currently set for placing emphasis predicted states over measured states
R = sig^2*eye(3); %Tune this, currently set for measurements being much more noisy than state predictions

% Need one extra slot because we store P_{k+1} after updating at step k
Parr = zeros(6,6,length(t)+1);
Parr(:,:,1) = cov0;

% Initial estimated state
xHat = zeros(6,length(t));
xHat(:,1) = [r0+ sig*randn(3,1); v0+ sig*randn(3,1)]; %Our initial estimate has some noise


%% ESTIMATION ----------------------------------------------------------

Estimate = EKF(t, Parr, xHat, y, H, Q, R, sig, cov0, muBennu, opts);
Estimate_BLS = EKF_BLS(t, Parr, xHat, y, H, Q, R, sig, cov0, muBennu, opts);

%NORMAL EKF METHOD
xHat = Estimate(1:6, :);
meas_hist = Estimate(7:9, :);
Pxx = Estimate(10,:);
Pyy = Estimate(11,:);
Pzz = Estimate(12,:);

%EKF + BLS METHOD
xHat_BLS = Estimate_BLS(1:6, :);
meas_hist_BLS = Estimate_BLS(7:9, :);
Pxx_BLS = Estimate_BLS(10,:);
Pyy_BLS = Estimate_BLS(11,:);
Pzz_BLS = Estimate_BLS(12,:);
%------------------------------------------------------------------------


%% Plotting -------------------------------------------------------------

% Precompute useful error quantities
err_EKF     = y(:,1:3)' - xHat(1:3,:);
err_BLS     = y(:,1:3)' - xHat_BLS(1:3,:);
err_meas    = y(:,1:3)' - meas_hist;
err_measBLS = y(:,1:3)' - meas_hist_BLS;

errNorm_EKF = vecnorm(err_EKF);
errNorm_BLS = vecnorm(err_BLS);

sig3_EKF = 3*[sqrt(Pxx); sqrt(Pyy); sqrt(Pzz)];
sig3_BLS = 3*[sqrt(Pxx_BLS); sqrt(Pyy_BLS); sqrt(Pzz_BLS)];

%% 1) 3D Trajectory Comparison
figure;
plot3(y(:,1), y(:,2), y(:,3), 'k-', 'LineWidth', 1.6); hold on;
plot3(xHat(1,:), xHat(2,:), xHat(3,:), 'b--', 'LineWidth', 1.3);
plot3(xHat_BLS(1,:), xHat_BLS(2,:), xHat_BLS(3,:), 'r-.', 'LineWidth', 1.3);

fv = stlread('g_06290mm_spc_obj_0000n00000_v008.stl');
patch('Faces', fv.ConnectivityList, ...
      'Vertices', fv.Points, ...
      'FaceColor', [0.8 0.8 0.8], ...
      'EdgeColor', 'none', ...
      'FaceLighting', 'gouraud');

axis equal
grid on
camlight headlight
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
title('Truth and Estimated Trajectories Around Bennu');
legend('Truth', 'EKF', 'EKF + BLS', 'Bennu Shape', 'Location', 'best');

%% 2) Position Estimation Error Components
figure;
labels = {'X Error (m)', 'Y Error (m)', 'Z Error (m)'};
titles = {'Position Error in X', 'Position Error in Y', 'Position Error in Z'};

for k = 1:3
    subplot(3,1,k)
    plot(t, err_EKF(k,:), 'b-', 'LineWidth', 1.3); hold on;
    plot(t, err_BLS(k,:), 'r--', 'LineWidth', 1.3);
    yline(0,'k:');
    grid on
    xlabel('Time (s)');
    ylabel(labels{k});
    title(titles{k});
    legend('EKF', 'EKF + BLS', 'Zero Error', 'Location', 'best');
end
sgtitle('Componentwise Position Estimation Errors');

%% 3) Measurement Error vs Estimate Error
figure;
for k = 1:3
    subplot(3,1,k)
    scatter(t, err_meas(k,:), 10, 'b', 'filled'); hold on;
    scatter(t, err_measBLS(k,:), 10, 'r', 'filled');
    plot(t, err_EKF(k,:), 'b-', 'LineWidth', 1.4);
    plot(t, err_BLS(k,:), 'r-', 'LineWidth', 1.4);
    yline(0,'k:');
    grid on
    xlabel('Time (s)');
    ylabel(labels{k});
    title(['Measurement and Estimation Error in ', char('X'+k-1)]);
    legend('Raw Meas Error (EKF)', 'BLS Meas Error', ...
           'Est Error (EKF)', 'Est Error (EKF+BLS)', ...
           'Zero Error', 'Location', 'best');
end
sgtitle('Measurement Errors vs Filtered Estimation Errors');

%% 4) Error Norm Comparison
figure;
plot(t, errNorm_EKF, 'b-', 'LineWidth', 1.6); hold on;
plot(t, errNorm_BLS, 'r--', 'LineWidth', 1.6);
grid on
xlabel('Time (s)');
ylabel('Position Error Norm (m)');
title('Norm of Position Estimation Error');
legend('EKF', 'EKF + BLS', 'Location', 'best');

%% 5) Improvement of EKF+BLS over EKF
figure;
plot(t, errNorm_EKF - errNorm_BLS, 'k', 'LineWidth', 1.5);
yline(0,'r--','LineWidth',1.2);
grid on
xlabel('Time (s)');
ylabel('\Delta Error Norm (m)');
title('Improvement from EKF+BLS Relative to EKF');
legend('EKF Error Norm - EKF+BLS Error Norm', 'Equal Performance', 'Location', 'best');

%% 6) True Error vs 3-Sigma Bounds (Consistency Check)
figure;
for k = 1:3
    subplot(3,1,k)
    plot(t, err_EKF(k,:), 'b-', 'LineWidth', 1.2); hold on;
    plot(t,  sig3_EKF(k,:), 'b--', 'LineWidth', 1.1);
    plot(t, -sig3_EKF(k,:), 'b--', 'LineWidth', 1.1);

    plot(t, err_BLS(k,:), 'r-', 'LineWidth', 1.2);
    plot(t,  sig3_BLS(k,:), 'r--', 'LineWidth', 1.1);
    plot(t, -sig3_BLS(k,:), 'r--', 'LineWidth', 1.1);

    yline(0,'k:');
    grid on
    xlabel('Time (s)');
    ylabel(labels{k});
    title(['Error and 3\sigma Bounds in ', char('X'+k-1)]);
    legend('EKF Error', 'EKF +3\sigma', 'EKF -3\sigma', ...
           'EKF+BLS Error', 'BLS +3\sigma', 'BLS -3\sigma', ...
           'Zero Error', 'Location', 'best');
end
sgtitle('Consistency Check: Actual Error vs Predicted 3\sigma Bounds');

%% 7) Covariance Comparison Only
figure;
subplot(3,1,1)
plot(t, Pxx, 'b-', 'LineWidth', 1.3); hold on;
plot(t, Pxx_BLS, 'r--', 'LineWidth', 1.3);
grid on
xlabel('Time (s)');
ylabel('Variance (m^2)');
title('X Error Covariance');
legend('EKF', 'EKF+BLS', 'Location', 'best');

subplot(3,1,2)
plot(t, Pyy, 'b-', 'LineWidth', 1.3); hold on;
plot(t, Pyy_BLS, 'r--', 'LineWidth', 1.3);
grid on
xlabel('Time (s)');
ylabel('Variance (m^2)');
title('Y Error Covariance');
legend('EKF', 'EKF+BLS', 'Location', 'best');

subplot(3,1,3)
plot(t, Pzz, 'b-', 'LineWidth', 1.3); hold on;
plot(t, Pzz_BLS, 'r--', 'LineWidth', 1.3);
grid on
xlabel('Time (s)');
ylabel('Variance (m^2)');
title('Z Error Covariance');
legend('EKF', 'EKF+BLS', 'Location', 'best');

sgtitle('Estimated Error Covariance Comparison');

%% 8) RMS Summary in Command Window
measErrNorm_EKF = sqrt(sum(err_meas.^2, 1, 'omitnan'));
measErrNorm_BLS = sqrt(sum(err_measBLS.^2, 1, 'omitnan'));

% Remove times where there was no measurement at all
measErrNorm_EKF = measErrNorm_EKF(~isnan(measErrNorm_EKF));
measErrNorm_BLS = measErrNorm_BLS(~isnan(measErrNorm_BLS));

rms_meas_EKF  = sqrt(mean(measErrNorm_EKF.^2));
rms_meas_BLS  = sqrt(mean(measErrNorm_BLS.^2));
rms_est_EKF   = sqrt(mean(errNorm_EKF(7:size(t,1)).^2)); %Removed the errors before our first measurement since it was skewing the error dramatically
rms_est_BLS   = sqrt(mean(errNorm_BLS(7:size(t,1)).^2)); %Removed the errors before our first measurement since it was skewing the error dramatically

fprintf('\n------ RMS Error Summary ------\n');
fprintf('Raw Measurement RMS Error (EKF):      %.4f m\n', rms_meas_EKF);
fprintf('BLS Measurement RMS Error:            %.4f m\n', rms_meas_BLS);
fprintf('EKF Estimate RMS Error:               %.4f m\n', rms_est_EKF);
fprintf('EKF + BLS Estimate RMS Error:         %.4f m\n', rms_est_BLS);
fprintf('RMS Improvement from BLS:             %.4f m\n', rms_est_EKF - rms_est_BLS);
fprintf('Percent Improvement from BLS:         %.2f %%\n', ...
    100*(rms_est_EKF - rms_est_BLS)/rms_est_EKF);

%% Functions---------------------------------------------------------------
%--------------------------------------------------------------------------
%Dynamics
function dx = bennuProp(t,x,muBody,aNoise)
    if nargin < 4
        aNoise = [0;0;0];
    end

    % Constants
    r = norm(x(1:3)); v = norm(x(4:6)); phi = reshape(x(7:42), [6 6]);
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

    a = -muBody * x(1:3) / r^3 + j2BCI + a_sun + a_SR + aNoise; %All perturbations

    A = zeros(6,6); A(1:3,4:6) = eye(3);
    A(4:6,1:3) = dfdr+ dj2dr + daSun_dr + daSRP_dr;  %State matrix with perturbations

    phiDot = A * phi; %Derivative of STM
    
    dx = [x(4:6); a; reshape(phiDot, [36 1])];

end

%Normal EKF Estimator
function estimate = EKF(t, Parr, xHat, y, H, Q, R, sig, cov0, muBennu, opts)
    for i = 1:length(t)-1
        % Current updated covariance/state
        P = Parr(:,:,i);
        xPlus = xHat(:,i);

        % --- Prediction step ---
        % Propagate estimated state over one time step, together with STM
        % This way our predictions are based off previous estimates
        xAug0 = [xPlus; reshape(eye(6),[36 1])];
        [~, yProp] = ode45(@(tt,xx) bennuProp(tt,xx,muBennu), [t(i) t(i+1)], xAug0, opts);

        xBar = yProp(end,1:6)';                     % predicted state
        phi  = reshape(yProp(end,7:42),[6 6]);     % predicted STM from estimated state
        PBar = phi*P*phi' + Q;                      % predicted error covariance

        % --- Measurement update ---
        if mod(i,6) == 0   % measurement every 10 min since dt = 100 s
            meas = y(i+1,1:3)' + sig*randn(3,1); %Measurements are true location + noise
            meas_hist(:,i+1) = meas;

            b = meas - H*xBar;
            K = PBar * H' / (H*PBar*H' + R);    %Kalman gains

            xHat(:,i+1) = xBar + K*b;           %Uses prediction and measurement
            Phat(:,:,i+1) = PBar - K*H*PBar;    %Uses prediction and measurement
        else
            meas_hist(:,i+1) = [NaN; NaN; NaN];   % no measurement this step

            xHat(:,i+1) = xBar;       % prediction only
            Phat(:,:,i+1) = PBar;     % covariance prediction only
        end

        Pxx(i+1) = Phat(1,1,i+1);
        Pyy(i+1) = Phat(2,2,i+1);
        Pzz(i+1) = Phat(3,3,i+1);

        % carry updated covariance forward
        Parr(:,:,i+1) = Phat(:,:,i+1);
    end
    % initialize first entries for plotting consistency
    meas_hist(:,1) = [NaN; NaN; NaN];
    Pxx(1) = cov0(1,1);
    Pyy(1) = cov0(2,2);
    Pzz(1) = cov0(3,3);

    estimate = [xHat; meas_hist; Pxx; Pyy; Pzz];
end

%EKF Modified by Batch Least Squares Estimator
function estimate = EKF_BLS(t, Parr, xHat, y, H, Q, R, sig, cov0, muBennu, opts)
    for i = 1:length(t)-1
        % Current updated covariance/state
        P = Parr(:,:,i);
        xPlus = xHat(:,i);

        % --- Prediction step ---
        % Propagate estimated state over one time step, together with STM
        % This way our predictions are based off previous estimates
        xAug0 = [xPlus; reshape(eye(6),[36 1])];
        [~, yProp] = ode45(@(tt,xx) bennuProp(tt,xx,muBennu), [t(i) t(i+1)], xAug0, opts);

        xBar = yProp(end,1:6)';                     % predicted state
        phi  = reshape(yProp(end,7:42),[6 6]);     % predicted STM from estimated state
        PBar = phi*P*phi' + Q;                      % predicted error covariance

        % --- Measurement update ---
        if mod(i,6) == 0   % update every 600 s since dt = 100 s
            meas_num = 100;

            meas_batch = zeros(3, meas_num);
            for j = 1:meas_num
                meas_batch(:,j) = y(i+1,1:3)' + sig*randn(3,1);
            end

            [x_bls_meas, R_bls] = BLS(meas_batch, R);
            meas_hist(:,i+1) = x_bls_meas;

            b = x_bls_meas - H*xBar;
            K = PBar * H' / (H*PBar*H' + R_bls);

            xHat(:,i+1) = xBar + K*b;
            Phat(:,:,i+1) = PBar - K*H*PBar;
        else
            meas_hist(:,i+1) = [NaN; NaN; NaN];

            xHat(:,i+1) = xBar;
            Phat(:,:,i+1) = PBar;
        end

        Pxx(i+1) = Phat(1,1,i+1);
        Pyy(i+1) = Phat(2,2,i+1);
        Pzz(i+1) = Phat(3,3,i+1);

        % carry updated covariance forward
        Parr(:,:,i+1) = Phat(:,:,i+1);
    end
    % initialize first entries for plotting consistency
    meas_hist(:,1) = [NaN; NaN; NaN];
    Pxx(1) = cov0(1,1);
    Pyy(1) = cov0(2,2);
    Pzz(1) = cov0(3,3);

    estimate = [xHat; meas_hist; Pxx; Pyy; Pzz];
end

%Batch Least Squares Estimator
function [estimate, R_bls] = BLS(meas_batch, R)
    m = size(meas_batch,2);

    Y = reshape(meas_batch, [3*m,1]);
    Hbig = kron(ones(m,1), eye(3));
    Rbig = kron(eye(m), R);

    estimate = (Hbig' / Rbig * Hbig) \ (Hbig' / Rbig * Y);
    R_bls = inv(Hbig' / Rbig * Hbig);
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
