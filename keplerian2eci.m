function [pos, vel] = keplerian2eci(a, e, i, raan, aop, ta, muValue)

    % Convert Keplerian orbital elements to Cartesian coordinates
    % All distance in km, all angles in radians

    % Constants
    %muValue = 3.986e5; % km^2 / s^2

    if nargin == 6
        muValue = 3.986e5; % km^2 / s^2
    end
    
    p = a * (1-e^2); % km
    h = sqrt(p*muValue);
    r = p / (1 + e * cos(ta)); % km

    %% Using 313 (raan, inc, aop+ta) - does not work, need to troubleshoot
    % RTN r and v
    % rRTN = r * [1;0;0];
    % vRTN = [e*sqrt(muEarth/p)*sin(ta); sqrt(muEarth/p)*(1+e*cos(ta)); 0];
    % rtn2eciDCM = orb2cart(raan, i, aop+ta);
    % pos = rtn2eciDCM * rRTN;
    % vel = rtn2eciDCM * vRTN;
    

    %% Using Perifocal frame: Keplerian -> Perifocal -> ECI
    % Perifocal r and v
    rPerifocal = r * [cos(ta); sin(ta); 0];
    vPerifocal = [-sqrt(muValue / p) * sin(ta);
                   sqrt(muValue / p) * (e + cos(ta));
                   0];
    % Rotate from perifocal to ECI
    r3aop = [cos(aop), -sin(aop), 0; sin(aop), cos(aop), 0; 0, 0, 1]; % Argument of periapsis
    r1inc = [1, 0, 0; 0, cos(i), -sin(i); 0, sin(i), cos(i)]; % Inclination
    r3raan = [cos(raan), -sin(raan), 0; sin(raan), cos(raan), 0; 0, 0, 1]; % RAAN

    R = r3raan * r1inc * r3aop; % DCM

    % Position and velocity in ECI
    pos = R * rPerifocal;
    vel = R * vPerifocal;
    
end