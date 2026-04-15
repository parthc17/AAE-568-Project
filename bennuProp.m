function dx = bennuProp(t,x,muBody,propSTM,aNoise)
    if nargin < 5
        aNoise = [0;0;0];
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

    a = -muBody * x(1:3) / r^3 + j2BCI + a_sun + a_SR + aNoise; %All perturbations

    if propSTM
        phi = reshape(x(7:42), [6 6]);
        A = zeros(6,6); A(1:3,4:6) = eye(3);
        A(4:6,1:3) = dfdr+ dj2dr + daSun_dr + daSRP_dr;  %State matrix with perturbations

        phiDot = A * phi; %Derivative of STM
    
        dx = [x(4:6); a; reshape(phiDot, [36 1])];
    else
        dx = [x(4:6); a];
    end

end