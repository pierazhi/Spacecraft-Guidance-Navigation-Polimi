% Spacecraft Guidance & Navigation
% Assignment # 2, Exercise 3
% Author: Piercarlo Fontana

%% 1) Check the visibility window 
clc; clear; close all;  set(0,'DefaultFigureVisible','on')
    
cspice_furnsh('/Users/pierazhi/Documents/MATLAB/GNC/Assignment 2 - Ex2/assignment02.tm')

% Define SGP4 settings
typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

rng('default');

r0 = [4307.844185282820, -1317.980749248651, 2109.210101634011];
v0 = [-0.110997301537882, -0.509392750828585, 0.815198807994189];

initial_epoch = '2024-11-18 16:30:00.000 UTC';
t0 = cspice_str2et(initial_epoch); 

final_epoch = '2024-11-18 20:30:00.000 UTC';
tf = cspice_str2et(final_epoch); 

mu = cspice_bodvrd('Moon', 'GM', 1);
r_Moon = cspice_bodvrd('Moon', 'RADII', 3);
r_Moon = r_Moon(1);

P0 = diag([10,1,1,0.001,0.001,0.001,0.00001,0.00001]);

lander = 'MOONLANDER';
LAT = deg2rad(78);
LON = deg2rad(15);
ALT = 0;

radii = cspice_bodvrd('MOON', 'RADII', 3);
re = radii(1); rp = radii(2);
flat = (re - rp) / re;
tspan = t0:30:tf;

[xf1, tf1, tt1, xx1] = keplerian_propagator(tspan, [r0 v0], 'Moon', 'keplerian');

%figure
% plot3(xx1(:, 1), xx1(:, 2), xx1(:, 3), 'b', 'LineWidth', 1.5) % Plot trajectory
% hold on

% Plot starting and arrival points
% scatter3(xx1(1, 1), xx1(1, 2), xx1(1, 3), 50, 'r', 'filled')  % Starting point
% scatter3(xx1(end, 1), xx1(end, 2), xx1(end, 3), 50, 'g', 'filled') % Arrival point
% [x, y, z] = sphere(20);
% surf(re*x, re*y, re*z, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]); % Gray Moon
% legend('Trajectory', 'Starting Point', 'Arrival Point', 'Moon', 'Location', 'best')
% grid on
% axis equal
% xlabel('X (km)')
% ylabel('Y (km)')
% zlabel('Z (km)')
% title('Orbiter orbit @MCMF')

pos_lander = cspice_pgrrec('MOON', LON, LAT, 0, re, flat);
fprintf('-- Lander Position --\n');
fprintf('X: %.3f\n', pos_lander(1))
fprintf('Y: %.3f\n', pos_lander(2))
fprintf('Z: %.3f\n\n', pos_lander(3))

for i = 1:length(tspan)
    [orbiter_range(i),orbiter_elevation(i)]  = lander_rae_first_point(LAT, LON, xx1(i, :), tspan(i));
end

% Convert ET to UTC and datenum for plotting
utc_dates = arrayfun(@(t) cspice_et2utc(t, 'ISOC', 3), tspan, 'UniformOutput', false);
date_nums = datenum(utc_dates, 'yyyy-mm-ddTHH:MM:SS.FFF');

% Find visibility indices (when elevation > 0)
visibility_indices = find(orbiter_elevation > deg2rad(0));

% Check if visibility indices are non-empty
if ~isempty(visibility_indices)
    % Print the start and end visibility times using curly braces to access the cell array elements
    fprintf('Visibility Period: Start %s - End %s\n', datestr(date_nums(visibility_indices(1)), 'dd-mm-yyyy HH:MM'), datestr(date_nums(visibility_indices(end))));
else
    fprintf('No visibility period found.\n');
end

figure;
scatter(date_nums, rad2deg(orbiter_elevation), 40, 'red', 'Marker', '.');
hold on;
title(['Elevation @', lander]);
xlabel('Date');
ylabel('Elevation [deg]');
ylim([rad2deg(min(orbiter_elevation)), rad2deg(max(orbiter_elevation)) + 2]); % Adjust ylim to show data near min elevation
legend('Elevation', 'Location', 'Best');
grid on;
datetick('x', 'dd-mm-yyyy HH:MM', 'keepticks'); % Format as readable dates


%Create the figure for range plots
% figure();
% scatter(date_nums, orbiter_range, 40, 'red', 'Marker', '.');
% hold on;
% title(['Range @', lander]);
% xlabel('Date');
% ylabel('Range [km]');
% grid on;
% datetick('x', 'dd-mm-yyyy HH:MM:SS', 'keepticks');


%% 2) Simulate Measurements
clc;

% Preallocate arrays
range = zeros(1, length(tspan));       % Range measurements
azimuth = zeros(1, length(tspan));     % Azimuth measurements
elevation = zeros(1, length(tspan));   % Elevation measurements

meas_noise = 0.1;

for i = 1:length(tspan)

    [range(i), azimuth(i), elevation(i)] = lander_rae(lander, tspan(i), xx1(i, :));
    range_noise(i) = mvnrnd(range(i), meas_noise^2);
    xx_noise(i, 1:3) = mvnrnd(xx1(i, 1:3), diag([meas_noise^2 meas_noise^2 meas_noise^2]));

end

% Convert ET to UTC and datenum for plotting
utc_dates2 = arrayfun(@(t) cspice_et2utc(t, 'ISOC', 3), tspan, 'UniformOutput', false);
date_nums = datenum(utc_dates2, 'yyyy-mm-ddTHH:MM:SS.FFF');

% Compute error
noisypos_error = xx1(:, 1:3) - xx_noise(:, 1:3); % Error in position components
noisyrange_error = range - range_noise; % Error in range

% Plot position error over time
figure;
plot(date_nums, noisypos_error(:, 1), '.'); hold on;
plot(date_nums, noisypos_error(:, 2), '.');
plot(date_nums, noisypos_error(:, 3), '.');
yline(3*meas_noise, 'r--'); yline(-3*meas_noise, 'r--');
title('Position Errors Over Time');
xlabel('Date');
ylabel('Error [km]');
legend('X Error', 'Y Error', 'Z Error', '±3σ');
grid on;
datetick('x', 'dd-mm-yyyy HH:MM', 'keepticks');

% Plot range error over time
figure;
plot(date_nums, noisyrange_error, '.'); hold on;
yline(3*meas_noise, 'r--'); yline(-3*meas_noise, 'r--');
title('Range Errors Over Time');
xlabel('Date');
ylabel('Error [km]');
legend('Error', '±3σ');
grid on;
datetick('x', 'dd-mm-yyyy HH:MM', 'keepticks');

% Plot evolution of position vector in time
figure;
hold on;
plot(date_nums, xx1(:, 1), 'r','LineWidth', 2);
plot(date_nums, xx1(:, 2), 'g','LineWidth', 2)
plot(date_nums, xx1(:, 3), 'b','LineWidth', 2)
title('Orbiter Position Components Time Evolution @ MCIF');
xlabel('Date');
ylabel('[km]');
legend('X', 'Y', 'Z', 'Location', 'Best');
grid on;
datetick('x', 'dd-mm-yyyy HH:MM', 'keepticks');
    
figure;
hold on;
plot(date_nums, range, 'r','LineWidth', 2);
title(['Range of ', lander]);
xlabel('Date');
ylabel('Range [km]');
legend('Range', 'Location', 'Best');
grid on;
datetick('x', 'dd-mm-yyyy HH:MM', 'keepticks');
   
%% 3) Estimate the lunar orbiter absolute state
clc; 

x0 = [r0 v0];
x0_noisy = mvnrnd(x0, P0(1:6, 1:6))';

n = 6;
alpha = 0.01;
beta = 2;

c = alpha^2 * n;
chi_0 = x0_noisy;
W0_m = 1 - n / (alpha^2 * n);
W0_c = (2 - alpha^2 + beta) - n / (alpha^2*n);
Wi = (1 / (2*alpha^2*n))*ones(1, 2*n);

W_m = [W0_m Wi];    
W_c = [W0_c Wi];

Rk = 0.1^2 * diag(ones(1, 3));

x_est(:, 1) = x0_noisy;
P_est(:, :, 1) = P0(1:6, 1:6);

for k = 2:length(tspan)

        % --- Sigma Point Generation ---

        sigma_points(:,1) = x_est(:, k-1);  % First sigma point is the mean state
        dX0 = sqrtm(c * (P_est(:, :, k-1)));

        for i = 1:n
                sigma_points(:, i+1) = x_est(:, k-1) + dX0(:, i);  % Positive sigma points
                sigma_points(:, i+n+1) = x_est(:, k-1) - dX0(:, i);  % Negative sigma points
        end

         % --- Sigma Point Propagation & Simulate the measurements  ---

        for i = 1:2*n+1

            % Propagate each sigma point using the dynamics model

            [sigma_points(:, i), ~, ~, ~] = keplerian_propagator([tspan(k-1), tspan(k)], sigma_points(:, i), 'MOON', 'keplerian');
        end        

        gamma_points = sigma_points(1:3, :);

       % --- Prediction Step ---

        % Predicted state mean
        x_minus = sum(W_m .* sigma_points, 2);
        P_minus = W_c.*(sigma_points-x_minus)*(sigma_points-x_minus)';

        % Predicted measurement mean

        y_minus = sum(W_m .* gamma_points, 2);
        Pee_k = W_c.*(gamma_points - y_minus)*(gamma_points - y_minus)' + Rk;
        Pxy_k = W_c.*(sigma_points - x_minus)*(gamma_points - y_minus)';

       % --- Update Step ---
        K_k = Pxy_k / Pee_k;  % Kalman Gain
        x_plus = x_minus + K_k * (xx_noise(k, :)' - y_minus);  % Updated state
        P_plus = P_minus - K_k * Pee_k * K_k';  % Updated covariance

        x_est(:, k)    = x_plus;
        P_est(:, :, k) = P_plus;
        
        position_3sigma(:,k) = 3 * sqrt(diag(P_est(1:3, 1:3, k)));
        velocity_3sigma(:,k) = 3 * sqrt(diag(P_est(4:6, 4:6, k)));

end

position_x_error = abs(x_est(1, :) - xx1(:, 1)');
position_y_error = abs(x_est(2, :) - xx1(:, 2)');
position_z_error = abs(x_est(3, :) - xx1(:, 3)');
velocity_x_error = abs(x_est(4, :) - xx1(:, 4)');
velocity_y_error = abs(x_est(5, :) - xx1(:, 5)');
velocity_z_error = abs(x_est(6, :) - xx1(:, 6)');

start = 1;

figure
hold on
plot(date_nums(start:end), position_x_error(start:end),'r')
plot(date_nums(start:end), position_y_error(start:end),'g')
plot(date_nums(start:end), position_z_error(start:end),'b')
plot(date_nums(start:end), position_3sigma(3, start:end),'k')
grid on
legend('X Component', 'Y Component', 'Z Component', '3sigma Boundary')
xlabel('Date')
ylabel('Position Error [km]')
title('Time Evolution of The Position Error')
datetick('x', 'dd-mm-yyyy HH:MM', 'keepticks');
%magnifyOnFigure;

figure
hold on
plot(date_nums(start:end), velocity_x_error(start:end),'r')
plot(date_nums(start:end), velocity_y_error(start:end),'g')
plot(date_nums(start:end), velocity_z_error(start:end),'b')
plot(date_nums(start:end), velocity_3sigma(1, start:end),'k')
grid on
legend('VX Component', 'VY Component', 'VZ Component', '3sigma Boundary')
xlabel('Date')
ylabel('Velocity Error [km/s]')
title('Time Evolution of The Velocity Error')
datetick('x', 'dd-mm-yyyy HH:MM', 'keepticks');
%magnifyOnFigure;

%% 4) Estimate lander and orbiter
clc; clear gamma_points sigma_points

x0 = [r0 v0 LAT LON];   % Initial state
x0_noisy = mvnrnd(x0, P0)';  % Add noise to the initial state

n = 8;  % State dimension
alpha = 0.01;
beta = 2;

c = alpha^2 * n;
W0_m = 1 - n / (alpha^2 * n);
W0_c = (2 - alpha^2 + beta) - n / (alpha^2*n);
Wi = (1 / (2*alpha^2*n))*ones(1, 2*n);

W_m = [W0_m Wi];
W_c = [W0_c Wi];

Rk = 0.1^2 * diag(ones(1, 4));

x_est = zeros(8, length(tspan));
P_est = zeros(8, 8, length(tspan));

x_est(:, 1) = x0_noisy;  % Initialize estimated state
P_est(:, :, 1) = P0;     % Initialize covariance matrix

for k = 2:length(tspan)
    % --- Sigma Point Generation ---
    sigma_points_m(:, 1) = x_est(:, k-1); % First sigma point is the mean state
    dX0 = sqrtm(c * P_est(:, :, k-1));
    for i = 1:n
        sigma_points_m(:, i+1) = x_est(:, k-1) + dX0(:, i); % Positive sigma points
        sigma_points_m(:, i+n+1) = x_est(:, k-1) - dX0(:, i); % Negative sigma points
    end

    % --- Sigma Point Propagation ---
    for i = 1:2*n+1
        % Propagate orbiter state
        [x_final(1:6, i), ~, ~, ~] = keplerian_propagator([tspan(k-1), tspan(k)], sigma_points_m(1:6, i), 'MOON', 'keplerian');
        r_orbiter = x_final(1:3, i);

        % Compute lander position from LAT/LON
        lat = sigma_points_m(7, i);
        lon = sigma_points_m(8, i);

        range_4 = lander_rae_first_point(lat, lon, r_orbiter', tspan(k));
        % Form measurement vector for each sigma point
        gamma_points(:, i) = [r_orbiter; range_4];
    end
    
    sigma_points_prop = [x_final; sigma_points_m(7, :); sigma_points_m(8, :)];
    
    % --- Prediction Step ---
    x_minus = sum(W_m .* sigma_points_prop, 2);
    P_minus = W_c .* (sigma_points_prop - x_minus) * (sigma_points_prop - x_minus)';

    y_minus = sum(W_m .* gamma_points, 2);
    Pee_k = W_c .* (gamma_points - y_minus) * (gamma_points - y_minus)' + Rk;
    Pxy_k = W_c .* (sigma_points_prop - x_minus) * (gamma_points - y_minus)';

    % --- Update Step ---
    K_k = Pxy_k / Pee_k; % Kalman Gain
    x_plus = x_minus + K_k * ([xx_noise(k, :)'; range_noise(k)] - y_minus); % Updated state
    P_plus = P_minus - K_k * Pee_k * K_k'; % Updated covariance

    x_est(:, k) = x_plus;
    P_est(:, :, k) = P_plus;
    
    position_3sigma(:,k) = 3 * sqrt(diag(P_est(1:3, 1:3, k)));
    velocity_3sigma(:,k) = 3 * sqrt(diag(P_est(4:6, 4:6, k)));
end
position_x_error = abs(x_est(1, :) - xx1(:, 1)');
position_y_error = abs(x_est(2, :) - xx1(:, 2)');
position_z_error = abs(x_est(3, :) - xx1(:, 3)');
velocity_x_error = abs(x_est(4, :) - xx1(:, 4)');
velocity_y_error = abs(x_est(5, :) - xx1(:, 5)');
velocity_z_error = abs(x_est(6, :) - xx1(:, 6)');


start = 1;

figure
hold on
plot(date_nums(start:end), position_x_error(start:end),'r')
plot(date_nums(start:end), position_y_error(start:end),'g')
plot(date_nums(start:end), position_z_error(start:end),'b')
plot(date_nums(start:end), position_3sigma(2, start:end),'k')
grid on
legend('X Component', 'Y Component', 'Z Component', '3sigma Boundary')
xlabel('Date')
ylabel('Position Error [km]')
title('Time Evolution of The Position Error')
datetick('x', 'dd-mm-yyyy HH:MM', 'keepticks');
%magnifyOnFigure;

figure
hold on
plot(date_nums(start:end), velocity_x_error(start:end),'r')
plot(date_nums(start:end), velocity_y_error(start:end),'g')
plot(date_nums(start:end), velocity_z_error(start:end),'b')
plot(date_nums(start:end), velocity_3sigma(3, start:end),'k')
grid on
legend('VX Component', 'VY Component', 'VZ Component', '3sigma Boundary')
xlabel('Date')
ylabel('Velocity Error [km/s]')
title('Time Evolution of The Velocity Error')
datetick('x', 'dd-mm-yyyy HH:MM', 'keepticks');
%magnifyOnFigure;

% Save Figures

% saveAllFiguresToPath('/Users/pierazhi/Documents/MATLAB/GNC/Assignment 2 - Ex3/Images')

%% Functions

function [range, azimuth, elevation, rll_lander_sat] = lander_rae(lander, et, xx)
%--------------------------------------------------------------------------
% Function: lander_rae_first_point
% Author: Piercarlo Fontana
% Date: 30 May 2025
%
% Computes the initial range and elevation of a satellite as seen from a 
% specified location on the lunar surface, given by latitude and longitude.
% This function is typically used to evaluate the observability of a 
% satellite from a lander or surface observer at the first time step.
%
% Inputs:
%   lat  - Latitude of the observer on the Moon [rad]
%   lon  - Longitude of the observer on the Moon [rad]
%   xx   - Satellite state history matrix [nx6], where each row is [x y z vx vy vz]
%   tt   - Time vector corresponding to each state [nx1], in seconds past J2000 TDB
%
% Outputs:
%   range     - Slant range between the observer and the satellite at the first time step [km]
%   elevation - Elevation angle of the satellite above the local horizon at the first time step [rad]
%--------------------------------------------------------------------------

% Compute lander position in ECI
rv_lander = cspice_spkezr(lander, et, 'J2000', 'NONE', 'MOON');

% Compute lander-satellite vector in ECI
rv_lander_sat = [xx(1, 1:3) xx(1, 4:6)]' - rv_lander;

% Define station name
topoFrame = [lander, '_TOPO'];

% Transformation from ECI to topocentric frame
ROT_ECI2TOPO = cspice_sxform('J2000', topoFrame, et);

% Convert lander-satellite into topocentric frame
rv_lander_sat_topo = ROT_ECI2TOPO*rv_lander_sat;

% Compute range, azimuth and elevation using cspice_xfmsta
rll_lander_sat = cspice_xfmsta(rv_lander_sat_topo,'RECTANGULAR','LATITUDINAL','MOON');

range   = rll_lander_sat(1);   % [km]
azimuth = rll_lander_sat(2);   % [rad]
elevation = rll_lander_sat(3); % [rad]

end

function [range, elevation] = lander_rae_first_point(lat, lon, xx, tt)
%--------------------------------------------------------------------------
% Function: lander_rae_first_point
% Author: Piercarlo Fontana
% Date: 30 May 2025
%
% Computes the slant range and elevation of a satellite at the first time step
% as seen from a fixed location on the lunar surface, given by latitude and 
% longitude. The computation is performed using topocentric transformation 
% and SPICE routines for Moon-centered geometry.
%
% Inputs:
%   lat  - Observer's latitude on the Moon [rad]
%   lon  - Observer's longitude on the Moon [rad]
%   xx   - Satellite state matrix [nx6]; first three columns are position [km]
%   tt   - Scalar ephemeris time [s] corresponding to the first state, TDB past J2000
%
% Outputs:
%   range     - Distance from the lander to the satellite [km]
%   elevation - Elevation angle of the satellite above the local horizon [rad]
%--------------------------------------------------------------------------

radii = cspice_bodvrd('MOON', 'RADII', 3);
re = radii(1); rp = radii(2);
flat = (re - rp) / re;

r_lander = cspice_pgrrec('MOON', lon, lat, 0, re, flat);

ROT_ECI2IAU = cspice_pxform('J2000','IAU_MOON', tt); % rotation matrix from J200 to IAU Moon reference frame

IAU2TOPO = cspice_eul2m(lat - pi, pi - lon, pi/2, 2, 1, 2);

r_orbiter_IAU = ROT_ECI2IAU*xx(:, 1:3)';

diff_IAU = r_orbiter_IAU - r_lander;

diff_TOPO = IAU2TOPO*diff_IAU;

% diff_TOPO is the right one!

[range, ~, elevation] = cspice_reclat(diff_TOPO);

end

function [dxdt] = keplerian_rhs(t, x, GM, perturbed)
    %   Evaluates the right-hand-side of a 2-body (keplerian) propagator
    %   
    %
    % Author
    %   Name: ALESSANDRO 
    %   Surname: MORSELLI
    %   Research group: DART
    %   Department: DAER
    %   University: Politecnico di Milano 
    %   Creation: 24/10/2021
    %   Contact: alessandro.morselli@polimi.it
    %   Copyright: (c) 2021 A. Morselli, Politecnico di Milano. 
    %                  All rights reserved.
    %
    %
    % Notes:
    %   This material was prepared to support the course 'Satellite Guidance
    %   and Navigation', AY 2021/2022.
    %
    %
    % Inputs:
    %   t   : [ 1, 1] epoch (unused)
    %   x   : [6, 1] cartesian state vector wrt Solar-System-Barycentre and
    %                 State Transition Matrix elements
    %   GM  : [ 1, 1] gravitational constant of the body
    %
    % Outputs:
    %   dxdt   : [6,1] RHS, newtonian gravitational acceleration only
    %
    
    J2 = 0.0010826269;
    Re = cspice_bodvrd('Earth','RADII', 3);
    Re = Re(1);

    % Initialize right-hand-side
    dxdt = zeros(6,1);
    if isequal(perturbed, 'no')
            % Extract positions
            rr = x(1:3);
        
            % Compute square distance and distance
            dist2 = dot(rr, rr);
            dist = sqrt(dist2);
        
            % Position detivative is object's velocity
            dxdt(1:3) = x(4:6);   
            % Compute the gravitational acceleration using Newton's law
            dxdt(4:6) = - GM * rr /(dist*dist2);
    
    elseif isequal(perturbed, 'yes')       
            % Initialize right-hand-side and acceleration
            dxdt = zeros(6,1);
        
            % Extract positions
            rr = x(1:3);
        
            rotm = cspice_pxform('J2000', 'ITRF93', t);
            rr_ECEF = rotm * rr; 
        
            % Compute square distance and distance
            dist2 = dot(rr_ECEF, rr_ECEF);
            dist = sqrt(dist2);
        
            % Position detivative is object's velocity
            dxdt(1:3) = x(4:6);   
        
            % Compute the J2 acceleration 
            J2 = 0.0010826269;
            Re = cspice_bodvrd('Earth','RADII', 3);
            Re = Re(1);
            aa_J2 = 1.5*GM*J2*rr_ECEF/(dist*dist2)*(Re/dist)^2.*(5*(rr_ECEF(3)/dist)^2 - [1; 1; 3]);
            aa_J2 = rotm'*aa_J2;
           
            % Compute the gravitational acceleration using Newton's law and add
            % compute the total acceleration
            dxdt(4:6) = - GM * rr /(dist*dist2) + aa_J2;
    end
end

function [xf, tf, tt, xx] = keplerian_propagator(tspan, x0, attractor, perturbed)
%--------------------------------------------------------------------------
% Function: keplerian_propagator
% Author: Piercarlo Fontana
% Date: 30 May 2025
%
% Propagates a spacecraft's state using Keplerian or J2-perturbed dynamics.
%
% This function integrates the equations of motion using `ode113` and supports
% pure two-body or J2-perturbed motion depending on the selected mode.
%
% Inputs:
%   tspan     - Time span vector [t0, tf]
%   x0        - Initial state vector [x, y, z, vx, vy, vz]
%   attractor - String name of central body or its GM value
%   perturbed - Propagation model: 'keplerian' or 'J2'
%
% Outputs:
%   xf        - Final state vector [x, y, z, vx, vy, vz]
%   tf        - Final time of propagation
%   tt        - Time vector (nx1)
%   xx        - State history matrix (nx6)
%--------------------------------------------------------------------------

% Initialize propagation data
if isfloat(attractor)
        GM = attractor;
else
        GM = cspice_bodvrd(attractor, 'GM', 1);
end

options = odeset('reltol', 1e-12, 'abstol', [1e-6*ones(1,3),1e-9*ones(1,3)]);

if isequal(perturbed, 'keplerian')

    [tt, xx] = ode113(@(t,x) keplerian_rhs(t, x, GM, 'no'),tspan, x0, options);

elseif isequal(perturbed, 'J2')

    [tt, xx] = ode113(@(t,x) keplerian_rhs(t, x, GM, 'yes'),tspan, x0, options);
end

    % Extract state vector and Final Time
    xf = xx(end,1:6)';
    tf = tt(end);

end

function saveAllFiguresToPath(userDefinedPath)
% SAVEALLFIGURESTOPATH Saves all open figures to a user-defined path.
%
%   saveAllFiguresToPath(userDefinedPath) saves all open MATLAB figures to the
%   specified userDefinedPath. The function generates filenames as Figure_1.png,
%   Figure_2.png, etc.

% Input Validation
if nargin == 0
    error('User-defined path is required.');
end

if ~ischar(userDefinedPath) || isempty(userDefinedPath)
    error('User-defined path must be a non-empty string.');
end

if ~exist(userDefinedPath, 'dir')
    mkdir(userDefinedPath);
    disp(['Created directory: ' userDefinedPath]);
end

% List all open figures
openFigures = findall(0, 'Type', 'figure');

if isempty(openFigures)
    disp('No open figures found.');
    return
end

% Save figures
counter = 1;
for i = 1:numel(openFigures)
    currentFigure = openFigures(i);
    
    % Generate filename and full file path
    fileName = sprintf('Figure_%d', counter);
    fullFilePath = fullfile(userDefinedPath, [fileName '.png']);
    
    % Save the figure
    try
        saveas(currentFigure, fullFilePath);
        disp(['Figure ' num2str(i) ' saved to: ' fullFilePath]);
    catch ME
        disp(['Error saving figure ' num2str(i) ': ' ME.message]);
    end
    
    counter = counter + 1;
end
end