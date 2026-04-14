function motion = cartesian(t, x, muBody)
    
    % Constants
    muEarth = muBody;
    % if (muBody == null)
    %     muEarth = 3.986e5; % km^2 / s^2
    % end
    radEarth = 6378.1; % km
    j2 = 1.0826e-3; % Earth oblateness factor
    j3 = -2.53244e-6; % Earth j3 factor 

    area2massRatio = 5.4e-6; % km^2/kg
    G0 = 1.02e14; % kg*km/s^2
    dSun = 149e6; % km
    muSun = 1.32712440018e11;
    mmEarth = sqrt(muSun/dSun^3); %2*pi/(365*86400); 
    
    r = sqrt(x(1)^2 + x(2)^2 + x(3)^2); % Position magnitude

    % J2 Acceleration 
    j2leadingTerm = (-3*muEarth*j2*radEarth^2) / (2*r^5);
    j2ECI = j2leadingTerm * [x(1)*(1-5*(x(3)^2/r^2)); x(2)*(1-5*(x(3)^2/r^2)); x(3)*(3-5*(x(3)^2/r^2))];

    % J3 Acceleration
    j3termXY = 5*muEarth*radEarth^3*j3 / (2*r^7);
    j3termZ = muEarth*radEarth^3*j3 / (2*r^5);
    j3ECI = [j3termXY * (7*x(3)^2 / r^2 - 3) * x(3)*x(1); j3termXY * (7*x(3)^2 / r^2 - 3) * x(3)*x(2); j3termZ * (3 - 30*x(3)^2/r^2 + 35*x(3)^4/r^4)];

    % SRP Acceleration
    taSun = mmEarth*t;
    rSun = dSun * [cos(taSun); sin(taSun); 0];
    sHat = rSun/dSun;

    srpLeadingTerm = area2massRatio * G0 / (dSun^2);
    srpECI = srpLeadingTerm * sHat;
    
    % Third Body Acceleration
    rr = [x(1); x(2); x(3)];
    thirdBodyECI = -muSun .* (((rr-rSun) ./ (norm(rr-rSun)^3)) + (rSun ./ (norm(rSun)^3)));

    % Total Perturbing Acceleration
    lol = [0 0 0];
    perturbECI = lol;%[j2ECI + j3ECI + srpECI + thirdBodyECI];

    dxdt1 = x(4);
    dxdt2 = x(5);
    dxdt3 = x(6);
    dxdt4 = (- muEarth * x(1) / r^3) + perturbECI(1);
    dxdt5 = (- muEarth * x(2) / r^3) + perturbECI(2);
    dxdt6 = (- muEarth * x(3) / r^3) + perturbECI(3);

    motion = [dxdt1; dxdt2; dxdt3; dxdt4; dxdt5; dxdt6];

end
