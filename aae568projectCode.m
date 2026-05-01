clc; clear; close all;
% AAE 568 Project Script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SETUP CONSTANTS---------------------------------------------------------

% Solver options (shared between parts)
tol = 1e-12;
opts = odeset('RelTol', tol, 'AbsTol', tol);

% Shared Bennu / spacecraft constants
muBennu = 5.2;                 % m^3/s^2
radBennu = 250;                % m
j2Bennu = 3.9257534110e-2;
j3Bennu = 1.4711698072e-2;
j4Bennu = 3.0760445246e-2;

% Bennu inertia / mass properties
ixxBennu = 1.8130e9;
iyyBennu = 1.8836e9;
izzBennu = 2.0334e9;           % kg km^2
massBennu = 7.7e10;

% Shared local-orbit initial condition around Bennu
bennuOrbit0 = struct;
bennuOrbit0.a = 1.5e3;
bennuOrbit0.e = 0.05;
bennuOrbit0.i = deg2rad(45);
bennuOrbit0.raan = 0;
bennuOrbit0.w = 0;
bennuOrbit0.f = 0;

[r0, v0] = keplerian2eci( ...
    bennuOrbit0.a, bennuOrbit0.e, bennuOrbit0.i, ...
    bennuOrbit0.raan, bennuOrbit0.w, bennuOrbit0.f, muBennu);

x0Bennu = [r0; v0];
bennuPeriod = 2*pi*sqrt(bennuOrbit0.a^3/muBennu);   % s

% Shared process-noise setting for truth propagation
sig_a_truth = 1e-7;           % m/s^2

%--------------------------------------------------------------------------
% Part 1: Earth to Bennu transfer
%--------------------------------------------------------------------------
muSun_km = 1.32712440018e11;   % km^3/s^2
kmPerAU = 1.496192602435979E+08;             % km/AU

% Canonical heliocentric units
lStar = kmPerAU;                              % km
tStar = sqrt(lStar^3 / muSun_km);

% Nondimensional heliocentric constants
muSun = muSun_km / (lStar^3 / tStar^2);
g0 = (9.80665 / 1000) / (lStar / tStar^2);    % nondim
isp = 4000 / tStar;                           % nondim

earthState0 = struct;
earthState0.a    = 1.496657326987069E+08 / lStar;
earthState0.e    = 1.704313732350883E-02;
earthState0.i    = deg2rad(6.198205899446798E-03);
earthState0.raan = deg2rad(1.799051298362160E+02);
earthState0.w    = deg2rad(2.816178320319530E+02);
earthState0.f    = deg2rad(2.669012326892521E+02);

bennuHelioState0 = struct;
bennuHelioState0.a    = 1.684375473676060E+08 / lStar;
bennuHelioState0.e    = 2.037604127989827E-01;
bennuHelioState0.i    = deg2rad(6.032225719717332E+00);
bennuHelioState0.raan = deg2rad(1.949096214163919E+00);
bennuHelioState0.w    = deg2rad(6.640855196853940E+01);
bennuHelioState0.f    = deg2rad(2.693898971751065E+02);

% Transfer time span (nondimensional heliocentric time)
t0Transfer = 0;
tfTransfer = 8;
tspanTransfer = linspace(t0Transfer, tfTransfer, 80);

%--------------------------------------------------------------------------
% Part 2: Estimation around Bennu
%--------------------------------------------------------------------------
cov0 = diag([10, 10, 10, 0.1, 0.1, 0.1]);

dtEst = 100;                                  % s
tEst = 0:dtEst:(21.55*86400);                 % s
bennuPeriod_days = bennuPeriod / 86400;

%--------------------------------------------------------------------------
% Part 3: LQR station keeping around Bennu
%--------------------------------------------------------------------------
tfStation = 10 * bennuPeriod;                 % s
ntStation = 5000;
tStation = linspace(0, tfStation, ntStation);
dtStation = tfStation / ntStation;

nx = 6;
nu = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1 - Obtain trajectory to Bennu

[rBennu0,vBennu0] = keplerian2eci( ...
    bennuHelioState0.a, bennuHelioState0.e, bennuHelioState0.i, ...
    bennuHelioState0.raan, bennuHelioState0.w, bennuHelioState0.f, muSun);

[tBennu, xBennuTraj] = ode45(@(t,x) cartesian(t,x,muSun), ...
    tspanTransfer, [rBennu0;vBennu0], opts);

bennuStateFinal = xBennuTraj(end,:);

[rEarth0,vEarth0] = keplerian2eci( ...
    earthState0.a, earthState0.e, earthState0.i, ...
    earthState0.raan, earthState0.w, earthState0.f, muSun);

[tEarth, xEarthTraj] = ode45(@(t,x) cartesian(t,x,muSun), ...
    tspanTransfer, [rEarth0;vEarth0], opts);

bvpOptions = bvpset('Stats','on','RelTol',1e-1);

m0Transfer = 20;
uMaxTransfer = 0.6;

Initial_Guess = load('568ProjEarthToBennuGuess_1MASS.mat');
sol = Initial_Guess.sol_mass;
rhoTransfer = [20 10 5 2 1 0.1];

for k = 1:length(rhoTransfer)
    rho = rhoTransfer(k);
    sol = bvp4c(@(t,x) BVP_ode_mass(t,x,rho,uMaxTransfer), ...
                 @(ya,yb) BVP_BC_mass(ya,yb,rEarth0,vEarth0,bennuStateFinal,m0Transfer), ...
                 sol, bvpOptions);
end
sol_transfer = sol;

% Plot result
figure;
plot3(sol.y(1,:), sol.y(2,:), sol.y(3,:), 'LineWidth', 1.5); hold on
plot3(rEarth0(1),rEarth0(2),rEarth0(3),'o','MarkerSize',7)
plot3(bennuStateFinal(1),bennuStateFinal(2),bennuStateFinal(3),'rx','MarkerSize',7)
plot3(xBennuTraj(:,1),xBennuTraj(:,2),xBennuTraj(:,3))
plot3(xEarthTraj(:,1),xEarthTraj(:,2),xEarthTraj(:,3))
plot3(0,0,0,'go','MarkerSize',15)
axis equal
grid on
legend('Mass-inclusive transfer','Initial Position','Final Bennu Position','Bennu Orbit','Earth Orbit','Sun')
title('Earth-to-Bennu Transfer with Mass Dynamics')

% Plot mass history just to check that its working
t_mass = sol.x;
y_mass = sol.y;

m_hist = y_mass(7,:); 
mf_transfer = m_hist(end);

% Plot mass consumed
m0Transfer = m_hist(1);
mass_used = m0Transfer - m_hist;

% Print quick diagnostics
fprintf('Initial mass: %.6f\n', m_hist(1));
fprintf('Final mass:   %.6f\n', m_hist(end));
fprintf('Mass used:    %.6f\n', m_hist(1) - m_hist(end));
fprintf('Minimum mass over trajectory: %.6f\n', min(m_hist));

%Calculate Optimal Input for this trajectory
uTransfer = zeros(1, size(sol.y,2));
for i = 1:size(sol.y, 2)
    lr = sol.y(8:10,i);
    lv = sol.y(11:13, i);
    lm = sol.y(14,i);

    % Numerical floors
    lvNorm = max(norm(lv,2), 1e-8);
    m = sol.y(7,i);
    mEff   = max(m, 1e-6);

    % Control
    uHatStar  = -lv / lvNorm;
    S         = 1 + lv' * uHatStar / mEff - lm / (isp*g0);
    gammaStar = 0.5 * uMaxTransfer * (1 + tanh(-S / rho));
    uTransfer(i) = norm(gammaStar * uHatStar);
end

%Plot Input
figure()
plot(t_mass, uTransfer)
title('Optimal Input History')
xlabel('Time (nondim)')
ylabel('Input')
grid()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2 - Station Keeping with EKF+BLS State Feedback

% Use shared setup values
r0Station = x0Bennu(1:3);
v0Station = x0Bennu(4:6);

%% Reference trajectory
[~, yRefAug] = ode45( ...
    @(t,x) bennuProp(t, x, muBennu, true, [0;0;0], [0;0;0], 'central', true), ...
    tStation, ...
    [r0Station; v0Station; reshape(eye(6), [nx^2 1]); reshape(zeros(nx, nu), [nx*nu 1])], ...
    opts);

xRef = yRefAug(:,1:6);

phiHist   = matrixify(yRefAug(:, (nx+1):(nx+nx^2)), 6, 6);
BkIntHist = matrixify(yRefAug(:, (nx+nx^2+1):(nx+nx^2+nx*nu)), nx, nu);
BkHist    = pagemtimes(phiHist, BkIntHist);

%% LQR gain computation
Q_lqr  = diag([1, 1, 1, 0.5, 0.5, 0.5]);
R_lqr  = diag([1, 1, 1]);
Qf_lqr = 2 * Q_lqr;

KHist = zeros(nu, nx, ntStation);
Pk = Qf_lqr;

for k = ntStation-1:-1:1
    Ak = phiHist(:,:,k+1) / phiHist(:,:,k);
    Bk = BkHist(:,:,k+1) - Ak * BkHist(:,:,k);

    S = R_lqr + Bk' * Pk * Bk;
    KHist(:,:,k) = S \ (Bk' * Pk * Ak);

    Pk = Q_lqr + Ak' * Pk * Ak - Ak' * Pk * Bk * (S \ (Bk' * Pk * Ak));
end

%% EKF + BLS settings
sigMeas = 2;
measStride = 2;
measBatchCount = 1:2:101;

H_meas = [eye(3) zeros(3,3)];
Q_ekf = 5*diag([1e-6, 1e-6, 1e-6, 1e-10, 1e-10, 1e-10]);
R_raw = sigMeas^2 * eye(3);

numCases = numel(measBatchCount);
maxBatch = max(measBatchCount);

% Fixed random realizations for fair comparison
rng(1);
processNoiseAll = sig_a_truth * randn(3, ntStation-1);
%processNoiseAll = kron(ones(1, ntStation-1), sig_a_truth * randn(3, 1));
measNoiseAll = sigMeas * randn(3, maxBatch, ntStation);

% Summary metrics
runtimeTotal = zeros(1, numCases);
runtimeBLS = zeros(1, numCases);
rmsEstCases = zeros(1, numCases);
rmsRawCases = zeros(1, numCases);
rmsBLSCases = zeros(1, numCases);

for iCase = 1:numCases
    batchSize = measBatchCount(iCase);

    % Reinitialize states for this case
    uRef = zeros(nu, ntStation);
    xTruth = zeros(nx, ntStation);
    xEst   = zeros(nx, ntStation);
    uHist  = zeros(nu, ntStation-1);

    rawMeasHist = NaN(3, ntStation);
    blsMeasHist = NaN(3, ntStation);
    PHist = zeros(nx, nx, ntStation);
    PHist(:,:,1) = cov0;

    xTruth(:,1) = [r0Station + [10;10;10]; v0Station];
    xEst(:,1)   = [r0Station + measNoiseAll(:,1,1); v0Station];

    blsTimeAccum = 0;
    totalTimer = tic;
    ODE_timer_accum = 0;
    for k = 1:ntStation-1
        % Controller uses estimated state
        errEst_k = xEst(:,k) - xRef(k,:)';
        u_k = uRef(:,k) - KHist(:,:,k) * errEst_k;
        uHist(:,k) = u_k;

        % Truth propagation with fixed pre-generated process noise
        ODE_timer = tic;
        aNoise_k = processNoiseAll(:,k);
        [~, yTruthSeg] = ode45( ...
            @(tt,xx) bennuProp(tt, xx, muBennu, false, aNoise_k, u_k, 'full', false), ...
            [tStation(k) tStation(k+1)], ...
            xTruth(:,k), ...
            opts);

        xTruth(:,k+1) = yTruthSeg(end,:)';

        % EKF prediction
        xAug0 = [xEst(:,k); reshape(eye(6), [36 1])];
        [~, yPredSeg] = ode45( ...
            @(tt,xx) bennuProp(tt, xx, muBennu, true, [0;0;0], u_k, 'full', false), ...
            [tStation(k) tStation(k+1)], ...
            xAug0, ...
            opts);
        ODE_timer_accum = ODE_timer_accum + toc(ODE_timer);
        xBar = yPredSeg(end,1:6)';
        phiPred = reshape(yPredSeg(end,7:42), [6 6]);
        PBar = phiPred * PHist(:,:,k) * phiPred' + Q_ekf;

        % Measurement update
        if mod(k, measStride) == 0
            measBatch = xTruth(1:3,k+1) + measNoiseAll(:,1:batchSize,k+1);

            rawMeasHist(:,k+1) = measBatch(:,1);

            blsTimer = tic;
            [zBLS, R_bls] = BLS(measBatch, R_raw);
            blsTimeAccum = blsTimeAccum + toc(blsTimer);

            blsMeasHist(:,k+1) = zBLS;

            innovation = zBLS - H_meas*xBar;
            KGain = PBar * H_meas' / (H_meas*PBar*H_meas' + R_bls);

            xEst(:,k+1) = xBar + KGain * innovation;
            PHist(:,:,k+1) = PBar - KGain * H_meas * PBar;
        else
            xEst(:,k+1) = xBar;
            PHist(:,:,k+1) = PBar;
        end
    end
    %ODE_timer_accum
    runtimeTotal(iCase) = toc(totalTimer);
    runtimeBLS(iCase) = blsTimeAccum;

    % Error metrics
    estErr = xTruth(1:3,:) - xEst(1:3,:);
    estErrNorm = vecnorm(estErr, 2, 1);

    rawMeasErr = xTruth(1:3,:) - rawMeasHist;
    blsMeasErr = xTruth(1:3,:) - blsMeasHist;

    validRaw = ~isnan(rawMeasErr(1,:));
    validBLS = ~isnan(blsMeasErr(1,:));

    rawNorm = vecnorm(rawMeasErr(:,validRaw), 2, 1);
    blsNorm = vecnorm(blsMeasErr(:,validBLS), 2, 1);

    rmsEstCases(iCase) = sqrt(mean(estErrNorm(validBLS).^2));
    rmsRawCases(iCase) = sqrt(mean(rawNorm.^2));
    rmsBLSCases(iCase) = sqrt(mean(blsNorm.^2));

    fprintf('Batch %3d | Total %.3f s | BLS %.3f s | RMS Est %.4f m\n', ...
        batchSize, runtimeTotal(iCase), runtimeBLS(iCase), rmsEstCases(iCase));
end

figure('Name','Batch Size Comparison');

subplot(3,1,1)
plot(measBatchCount, rmsEstCases, '-o', 'LineWidth', 1.5)
grid on
xlabel('measBatchCount')
ylabel('RMS Est Error (m)')
title('Estimation Accuracy vs Batch Size')

subplot(3,1,2)
plot(measBatchCount, runtimeBLS, '-o', 'LineWidth', 1.5)
grid on
xlabel('measBatchCount')
ylabel('BLS-only Runtime (s)') 
title('BLS Runtime vs Batch Size')

subplot(3,1,3)
plot(measBatchCount, runtimeTotal, '-o', 'LineWidth', 1.5)
grid on
xlabel('measBatchCount')
ylabel('Total Runtime (s)')
title('Total Runtime vs Batch Size')

figure('Name','Measurement Quality vs Batch Size');
plot(measBatchCount, rmsRawCases, '-o', 'LineWidth', 1.5, 'DisplayName', 'Raw RMS'); hold on
plot(measBatchCount, rmsBLSCases, '-s', 'LineWidth', 1.5, 'DisplayName', 'BLS RMS')
grid on
xlabel('measBatchCount')
ylabel('RMS Error (m)')
title('Raw vs BLS Measurement Error')
legend('Location','best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Open-loop truth for comparison
[~, yOpen] = ode45( ...
    @(t,x) bennuProp(t, x, muBennu, false, [0;0;0], [0;0;0], 'full', false), ...
    tStation, ...
    [r0Station + [10;10;10]; v0Station], ...
    opts);

%% Build outputs for plotting
yCtrl = xTruth.';
yEst  = xEst.';

trackErrTruth = yCtrl(:,1:6) - yRefAug(:,1:6);
trackErrEst   = yEst(:,1:6)  - yRefAug(:,1:6);

estErr = yCtrl(:,1:6)' - xEst;
estErrNorm = vecnorm(estErr(1:3,:), 2, 1);

sig3Hist = 3 * [ ...
    sqrt(squeeze(PHist(1,1,:)))'; ...
    sqrt(squeeze(PHist(2,2,:)))'; ...
    sqrt(squeeze(PHist(3,3,:)))' ];

%% Derived quantities for plotting
trackErrNorm = vecnorm(trackErrTruth', 2, 1);

posTruthErrNorm = vecnorm(trackErrTruth(:,1:3)', 2, 1);
posEstErrNorm   = vecnorm(estErr(1:3,:), 2, 1);

rawMeasErr = yCtrl(:,1:3)' - rawMeasHist;
blsMeasErr = yCtrl(:,1:3)' - blsMeasHist;

validRaw = ~isnan(rawMeasErr(1,:));
validBLS = ~isnan(blsMeasErr(1,:));

rawMeasErrNorm = NaN(1, ntStation);
blsMeasErrNorm = NaN(1, ntStation);
rawMeasErrNorm(validRaw) = vecnorm(rawMeasErr(:,validRaw), 2, 1);
blsMeasErrNorm(validBLS) = vecnorm(blsMeasErr(:,validBLS), 2, 1);

%% 1) 3D trajectory comparison (full view)
figure('Name','Trajectories: Full View');
plot3(yRefAug(:,1), yRefAug(:,2), yRefAug(:,3), 'LineWidth', 1.5, ...
    'DisplayName', 'Reference Orbit');
hold on
% plot3(yOpen(:,1), yOpen(:,2), yOpen(:,3), 'LineWidth', 1.5, ...
%     'DisplayName', 'Open-Loop Perturbed');
plot3(yCtrl(:,1), yCtrl(:,2), yCtrl(:,3), '--', 'LineWidth', 1.5, ...
    'DisplayName', 'Closed-Loop Truth');
plot3(yEst(:,1), yEst(:,2), yEst(:,3), ':', 'LineWidth', 1.5, ...
    'DisplayName', 'Closed-Loop Estimate');

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

%% 2) Tracking error components
figure('Name','Tracking Error Components');

subplot(2,1,1)
plot(tStation, trackErrTruth(:,1:3), 'LineWidth', 1.2)
grid on
xlabel('Time (s)')
ylabel('Position Error (m)')
title('Truth Position Tracking Error')
legend('e_x', 'e_y', 'e_z', 'Location', 'best')

subplot(2,1,2)
plot(tStation, trackErrTruth(:,4:6), 'LineWidth', 1.2)
grid on
xlabel('Time (s)')
ylabel('Velocity Error (m/s)')
title('Truth Velocity Tracking Error')
legend('e_{v_x}', 'e_{v_y}', 'e_{v_z}', 'Location', 'best')

%% 3) Norm summary plot
figure('Name','Norm Summaries');
plot(tStation, trackErrNorm, 'LineWidth', 1.4, 'DisplayName', 'Tracking Error Norm'); hold on
plot(tStation, posEstErrNorm, 'LineWidth', 1.4, 'DisplayName', 'Estimation Error Norm');
grid on
xlabel('Time (s)')
ylabel('Norm')
title('Tracking and Estimation Error Norms')
legend('Location', 'best')

%% 4) Control history
figure('Name','LQR Control History');
plot(tStation(1:end-1), uHist', 'LineWidth', 1.2)
grid on
xlabel('Time (s)')
ylabel('Control Acceleration (m/s^2)')
title('LQR Control Correction History')
legend('u_x', 'u_y', 'u_z', 'Location', 'best')

%% 5) Estimation error vs 3-sigma
figure('Name','Estimation Error vs 3-Sigma');

labels = {'X Error (m)', 'Y Error (m)', 'Z Error (m)'};
for idx = 1:3
    subplot(3,1,idx)
    plot(tStation, estErr(idx,:), 'LineWidth', 1.2); hold on
    plot(tStation,  sig3Hist(idx,:), '--', 'LineWidth', 1.0)
    plot(tStation, -sig3Hist(idx,:), '--', 'LineWidth', 1.0)
    yline(0, 'k:')
    grid on
    xlabel('Time (s)')
    ylabel(labels{idx})
    title(['Estimation Error and 3\sigma Bounds: ', char('X'+idx-1)])
    legend('Estimation Error', '+3\sigma', '-3\sigma', 'Zero Error', 'Location', 'best')
end

%% 6) Measurement quality: raw vs BLS norm
figure('Name','Measurement Error Norm');
plot(tStation, rawMeasErrNorm, '.', 'DisplayName', 'Raw Measurement Error Norm'); hold on
plot(tStation, blsMeasErrNorm, '.', 'DisplayName', 'BLS Measurement Error Norm');
grid on
xlabel('Time (s)')
ylabel('Position Error Norm (m)')
title('Raw vs BLS Measurement Error Norm')
legend('Location', 'best')

%% RMS summary
validRaw = ~isnan(rawMeasErr(1,:));
validBLS = ~isnan(blsMeasErr(1,:));

rawMeasErrNorm = vecnorm(rawMeasErr(:,validRaw), 2, 1);
blsMeasErrNorm = vecnorm(blsMeasErr(:,validBLS), 2, 1);

rmsRawMeas = sqrt(mean(rawMeasErrNorm.^2));
rmsBLSMeas = sqrt(mean(blsMeasErrNorm.^2));
rmsEst = sqrt(mean(estErrNorm(validBLS).^2));

fprintf('\n------ Station-Keeping + EKF/BLS Summary ------\n');
fprintf('Raw Measurement RMS Error:         %.4f m\n', rmsRawMeas);
fprintf('BLS Measurement RMS Error:         %.4f m\n', rmsBLSMeas);
fprintf('Closed-Loop Estimate RMS Error:    %.4f m\n', rmsEst);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3: Return Trip
tfStation = tfStation/tStar;
traveltime=14.5;
tfReturn = tfStation + traveltime;
tspanStation = linspace(0, tfStation, 80);
tspanReturn = linspace(tfStation, tfReturn, 80);
tspanReturn_Earth = linspace(0, tfReturn, 80);

[tBennu_R, xBennuTraj_R] = ode45(@(t,x) cartesian(t,x,muSun), ...
    tspanStation, [bennuStateFinal], opts);
[tEarth_R, xEarthTraj_R] = ode45(@(t,x) cartesian(t,x,muSun), ...
    tspanReturn_Earth, [rEarth0;vEarth0], opts);

bennuState0_R = xBennuTraj_R(end,:);
earthStateFinal_R = xEarthTraj_R(end,:);

Initial_Guess = load('Return_GuessFinal_1405days_umax061.mat');
sol = Initial_Guess.sol;
uMaxReturn = 0.61;
m0Return = mf_transfer;
for rho = 40:-10:10
    sol = bvp4c(@(t,x) BVP_ode_return_mass(t,x,rho,uMaxReturn), ...
                            @(ya,yb) BVP_BC_return_mass(ya,yb,bennuState0_R,earthStateFinal_R,m0Return), ...
                            sol, bvpOptions);
end
for rho = 10:-2:1
    sol = bvp4c(@(t,x) BVP_ode_return_mass(t,x,rho,uMaxReturn), ...
                            @(ya,yb) BVP_BC_return_mass(ya,yb,bennuState0_R,earthStateFinal_R,m0Return), ...
                            sol, bvpOptions);
end
for rho = 1:-.2:.1
    sol = bvp4c(@(t,x) BVP_ode_return_mass(t,x,rho,uMaxReturn), ...
                            @(ya,yb) BVP_BC_return_mass(ya,yb,bennuState0_R,earthStateFinal_R,m0Return), ...
                            sol, bvpOptions);
end
for rho = .1:-.01:.01
    sol = bvp4c(@(t,x) BVP_ode_return_mass(t,x,rho,uMaxReturn), ...
                            @(ya,yb) BVP_BC_return_mass(ya,yb,bennuState0_R,earthStateFinal_R,m0Return), ...
                            sol, bvpOptions);
end
sol_Return = sol;

% Plot result
figure;
hold on
plot3(rEarth0(1),rEarth0(2),rEarth0(3),'o','MarkerSize',7)
plot3(earthStateFinal_R(1),earthStateFinal_R(2),earthStateFinal_R(3),'go','MarkerSize',7)
plot3(bennuStateFinal(1),bennuStateFinal(2),bennuStateFinal(3),'rx','MarkerSize',7)
plot3(bennuState0_R(1),bennuState0_R(2),bennuState0_R(3),'bx','MarkerSize',7)
plot3(xBennuTraj(:,1),xBennuTraj(:,2),xBennuTraj(:,3))
plot3(xEarthTraj(:,1),xEarthTraj(:,2),xEarthTraj(:,3))
plot3(0,0,0,'go','MarkerSize',15)
axis equal
grid on
title('Earth-to-Bennu Transfer with Mass Dynamics')
% plot3(solinit.y(1,:), solinit.y(2,:), solinit.y(3,:))
plot3(sol.y(1,:), sol.y(2,:), sol.y(3,:))
windowsize = 2;
xlim([-windowsize,windowsize])
ylim([-windowsize, windowsize])
zlim([-windowsize,windowsize])
legend('Earth Pos (t=0)', 'Earth Pos at End of Return', 'Bennu Pos when we reach', 'Bennu Pos when we leave','Bennu Traj', 'Earth Traj' , 'Sun', 'Solved Trajectory')

uReturn = zeros(1, size(sol.y,2));
for i = 1:size(sol.y, 2)
    lr = sol.y(8:10, i);
    lv = sol.y(11:13, i);
    lm = sol.y(14, i);
    m  = sol.y(7, i);
    mEff = max(m, 1e-6);

    % Control input setup
    lvNorm = max(norm(lv,2), 1e-8);
    uHatStar  = -lv / lvNorm;
    S         = 1 + lv' * uHatStar / mEff - lm / (isp*g0);
    gammaStar = 0.5 * uMaxReturn * (1 + tanh(-S / rho));
    uReturn(i) = norm(gammaStar * uHatStar);
end

%Plot Input
figure()
plot(sol.x, uReturn)
title('Optimal Input History')
xlabel('Time (nondim)')
ylabel('Input')
grid()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function dx = bennuProp(t, x, muBody, propSTM, aNoise, u, modelType, propBk)
    % Consolidated Bennu propagator
    %
    % Inputs
    %   t         - time
    %   x         - state vector
    %   muBody    - Bennu gravitational parameter
    %   propSTM   - true/false, propagate STM
    %   aNoise    - 3x1 process-noise acceleration
    %   u         - 3x1 control acceleration
    %   modelType - 'central', 'j2', or 'full'
    %   propBk    - true/false, propagate Bk integral states (requires STM)
    %
    % State layouts supported:
    %   6-state:   [r; v]
    %   42-state:  [r; v; vec(Phi)]
    %   60-state:  [r; v; vec(Phi); vec(BkInt)]

    if nargin < 4 || isempty(propSTM)
        propSTM = false;
    end
    if nargin < 5 || isempty(aNoise)
        aNoise = [0;0;0];
    end
    if nargin < 6 || isempty(u)
        u = [0;0;0];
    end
    if nargin < 7 || isempty(modelType)
        modelType = 'full';
    end
    if nargin < 8 || isempty(propBk)
        propBk = false;
    end

    % Basic state
    rVec = x(1:3);
    vVec = x(4:6);
    r = norm(rVec);

    % Constants
    radBennu = 250;                 % m
    j2 = 3.9257534110e-2;

    % Base 2-body terms
    a = -muBody * rVec / r^3;
    dfdr = muBody * (3*rVec*rVec' / r^5 - eye(3)/r^3);

    % Start state Jacobian
    A = zeros(6,6);
    A(1:3,4:6) = eye(3);
    A(4:6,1:3) = dfdr;

    % -------------------------
    % J2 terms
    % -------------------------
    j2leadingTerm = (-3*muBody*j2*radBennu^2) / (2*r^5);
    j2BCI = j2leadingTerm * [ ...
        rVec(1)*(1-5*(rVec(3)^2/r^2));
        rVec(2)*(1-5*(rVec(3)^2/r^2));
        rVec(3)*(3-5*(rVec(3)^2/r^2)) ];

    dj2dx = -3*muBody*j2*radBennu^2 / (2*r^5) * ...
        [5*rVec(1)*(7*rVec(3)^2/r^2 - 1)/r^2 * rVec' + ...
         (1-5*rVec(3)^2/r^2) * [1,0,0] - ...
         10*rVec(1)*rVec(3)/r^2 * [0,0,1]];

    dj2dy = -3*muBody*j2*radBennu^2 / (2*r^5) * ...
        [5*rVec(2)*(7*rVec(3)^2/r^2 - 1)/r^2 * rVec' + ...
         (1-5*rVec(3)^2/r^2) * [0,1,0] - ...
         10*rVec(2)*rVec(3)/r^2 * [0,0,1]];

    dj2dz = -3*muBody*j2*radBennu^2 / (2*r^5) * ...
        [5*rVec(3)*(7*rVec(3)^2/r^2 - 3)/r^2 * rVec' + ...
         3*(1-5*rVec(3)^2/r^2) * [0,0,1]];

    dj2dr = [dj2dx; dj2dy; dj2dz];

    % -------------------------
    % Sun + SRP terms
    % -------------------------
    G = 6.674e-11;
    mSun = 1.989e30;
    mInAU = 1.496192602435979E+11;

    r_ben2sun = [1.12 * mInAU; 0; 0];
    d = r_ben2sun - rVec;

    a_sun = G*mSun * ( d / norm(d)^3 - r_ben2sun / norm(r_ben2sun)^3 );
    daSun_dr = G*mSun * ( 3*(d*d')/norm(d)^5 - eye(3)/norm(d)^3 );

    P_SR = 4.51*10^-6;
    c_R = 0.6;
    mSC = 20;
    A_exposed = 1;

    a_SRP = -P_SR*c_R*A_exposed*d/(mSC*norm(d));
    daSRP_dr = -P_SR*c_R*A_exposed/mSC * (d*d'/norm(d)^3 - eye(3)/norm(d));

    % -------------------------
    % Select model
    % -------------------------
    switch lower(modelType)
        case 'central'
            % already set from 2-body only

        case 'j2'
            a = a + j2BCI;
            A(4:6,1:3) = dfdr + dj2dr;

        case 'full'
            a = a + j2BCI + a_sun + a_SRP;
            A(4:6,1:3) = dfdr + dj2dr + daSun_dr + daSRP_dr;

        otherwise
            error('Unknown modelType: %s. Use ''central'', ''j2'', or ''full''.', modelType);
    end

    % Add noise and control
    a = a + aNoise + u;

    % -------------------------
    % Output by state size / flags
    % -------------------------
    if ~propSTM
        dx = [vVec; a];
        return
    end

    % STM case
    phi = reshape(x(7:42), [6 6]);
    phiDot = A * phi;

    if ~propBk
        dx = [vVec; a; reshape(phiDot, [36 1])];
        return
    end

    % STM + Bk integral case
    B = [zeros(3,3); eye(3)];
    BkDot = phi \ B;

    dx = [vVec;
          a;
          reshape(phiDot, [36 1]);
          reshape(BkDot, [18 1])];
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

function motion = BVP_ode_return_mass(t, x, rho, uMax)

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

function psi = BVP_BC_return_mass(ya, yb, bennuState0_R, earthStateFinal_R, m0Return)

    psi = [ya(1:6) - bennuState0_R(1:6)';   % initial position/velocity
           ya(7)   - m0Return;              % initial mass = final transfer mass
           yb(1:6) - earthStateFinal_R(1:6)'; % terminal position/velocity
           yb(14)];                         % free final mass => lambda_m(tf)=0
end
