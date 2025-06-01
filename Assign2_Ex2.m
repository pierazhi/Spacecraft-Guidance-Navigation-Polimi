% Spacecraft Guidance & Navigation
% Assignment # 2, Exercise 2
% Author: Piercarlo Fontana

%% 1) Compute visibility windows
clc; clearvars; cspice_kclear; set(0,'DefaultFigureVisible','on'); close all

cspice_furnsh('kernels/assignment02.tm')

% Define SGP4 settings
typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

longstr1 = '1 36036U 09059A   24323.76060260  .00000600  00000-0  20543-3 0  9995';
longstr2 = '2 36036  98.4396 148.4689 0001262 95.1025 265.0307 14.39727995790658';
    
% Constant for arcseconds to radians conversions
arcsec2rad = pi / (180*3600);
% Set gravitational parameters
mu = cspice_bodvrd('Earth', 'GM', 1);

rng('default');  % Set random seed for reproducibility

% Initialize the satrec structure, using the function twoline2rv
satrec = twoline2rv(longstr1, longstr2, typerun,'u', opsmode, whichconst);

% Get TLE epoch
[year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
et_ref = cspice_str2et(sat_epoch_str);

% Centuries from TDT 2000 January 1 00:00:00.000 
ttt = cspice_unitim(et_ref, 'ET', 'TDT')/cspice_jyear()/100;
ddpsi = -0.115178*arcsec2rad; %  [rad]
ddeps = -0.007395*arcsec2rad; %  [rad]

% Evaluate the TLE at the initial reference epoch
[satrec,rteme,vteme] = sgp4(satrec, 0.0);

% Get the osculating orbital elements
elts = cspice_oscelt( [rteme;vteme], et_ref, satrec.mu );

% Transform TEME to ECI vectors at reference epoch
ateme = [0;0;0];
[reci, veci, aeci] = teme2eci(rteme, vteme, ateme, ttt, ddpsi, ddeps);

% Start of the visibility window

sat_epoch_str_ini = '2024-11-18 20:30:00.000 UTC';  % directly in ISO 8601 format
et_0 = cspice_str2et(sat_epoch_str_ini);

% End of the visibility window

sat_epoch_str_fin = '2024-11-18 22:15:00.000 UTC';  % directly in ISO 8601 format
et_f = cspice_str2et(sat_epoch_str_fin);

npoints = round((et_0-et_ref)/60.0)+1;
tspan1 = linspace(et_ref, et_0, npoints);

% Propagate from TLE reference epoch to et_0 using sgp4

for i = 1:length(tspan1)

    % SGP4 propagation
    tsince = (tspan1(i) - et_ref)/60.0; % minutes from TLE epoch
    [~,rteme1,vteme1] = sgp4(satrec,  tsince);
    
    % Compute centuries from 2000-01-01T00:00:00.00 TDT
    ttt1 = cspice_unitim(tspan1(i), 'ET', 'TDT')/cspice_jyear()/100;
    
    % TEME to ECI conversion
    [reci1(:,i), veci1(:,i)] = teme2eci(rteme1, vteme1, [0.0;0.0;0.0], ttt1, ddpsi, ddeps);

end
reci1 = reci1';
veci1 = veci1';

% Equally spaced vector, with one minute time-step
npoints_60s = round((et_f-et_0)/60.0)+1;
tspan_visibility_60s = linspace(et_0, et_f, npoints_60s);
three.tspan{1} = tspan_visibility_60s;

npoints_30s = round((et_f-et_0)/30.0)+1;
tspan_visibility_30s = linspace(et_0, et_f, npoints_30s);
three.tspan{2} = tspan_visibility_30s;
three.tspan{3} = tspan_visibility_60s;

% Propagate from et_0 to et_f

for i = 1:length(tspan_visibility_30s)

    [~,rteme_30s,vteme_30s] = sgp4(satrec, (tspan_visibility_30s(i) - et_ref)/60.0);
    
    % Compute centuries from 2000-01-01T00:00:00.00 TDT
    ttt_30s = cspice_unitim(tspan_visibility_30s(i), 'ET', 'TDT')/cspice_jyear()/100;
    
    % TEME to ECI conversion
    [reci_30s(:,i), veci_30s(:,i)] = teme2eci(rteme_30s, vteme_30s, [0.0;0.0;0.0], ttt_30s, ddpsi, ddeps);

end
reci_30s = reci_30s';
veci_30s = veci_30s';

for i = 1:length(tspan_visibility_60s)

    [~,rteme_60s,vteme_60s] = sgp4(satrec,  (tspan_visibility_60s(i) - et_ref)/60.0);
    
    % Compute centuries from 2000-01-01T00:00:00.00 TDT
    ttt_60s = cspice_unitim(tspan_visibility_60s(i), 'ET', 'TDT')/cspice_jyear()/100;
    
    % TEME to ECI conversion
    [reci_60s(:,i), veci_60s(:,i)] = teme2eci(rteme_60s, vteme_60s, [0.0;0.0;0.0], ttt_60s, ddpsi, ddeps);

end
reci_60s = reci_60s';
veci_60s = veci_60s';

% figure
% plot3(reci_60s(:, 1), reci_60s(:, 2), reci_60s(:, 3))
% hold on
% scatter3(reci_60s(1, 1), reci_60s(1, 2), reci_60s(1, 3), 50, 'r', 'filled')
% scatter3(reci_60s(end, 1), reci_60s(end, 2), reci_60s(end, 3), 50, 'b', 'filled')
% legend('Trajectory', 'Starting Point', 'Arrival Point', 'Location', 'best')
% grid on
% title('SGP4 60s')
% 
% figure
% plot3(reci_30s(:, 1), reci_30s(:, 2), reci_30s(:, 3))
% hold on
% scatter3(reci_30s(1, 1), reci_30s(1, 2), reci_30s(1, 3), 50, 'r', 'filled')
% legend('Trajectory', 'Starting Point', 'Arrival Point', 'Location', 'best')
% grid on
% title('SGP4 30s')

% Keplerian Propagation

[xf1, tf1, tt1, xx1] = keplerian_propagator(tspan1, [reci veci], 'Earth', 'keplerian');
[xf_60s, ~, tt_60s, xx_60s] = keplerian_propagator(tspan_visibility_60s, xf1, 'Earth', 'keplerian');
[xf_30s, ~, tt_30s, xx_30s] = keplerian_propagator(tspan_visibility_30s, xf1, 'Earth', 'keplerian');

% figure
% plot3(xx_60s(:, 1), xx_60s(:, 2), xx_60s(:, 3))
% hold on
% scatter3(xx_60s(1, 1), xx_60s(1, 2), xx_60s(1, 3), 50, 'r', 'filled')
% scatter3(xx_60s(end, 1), xx_60s(end, 2), xx_60s(end, 3), 50, 'b', 'filled')
% legend('Trajectory', 'Starting Point', 'Arrival Point', 'Location', 'best')
% grid on
% title('Keplerian 60s')
% 
% figure
% plot3(xx_30s(:, 1), xx_30s(:, 2), xx_30s(:, 3))
% hold on
% scatter3(xx_30s(1, 1), xx_30s(1, 2), xx_30s(1, 3), 50, 'r', 'filled')
% legend('Trajectory', 'Starting Point', 'Arrival Point', 'Location', 'best')
% grid on
% title('Keplerian 30s')

diff_x_60s = norm(xx_60s(:, 1) - reci_60s(:, 1));
diff_x_30s = norm(xx_30s(:, 1) - reci_30s(:, 1));   

%% 2 & 3) Simulate Measurements Over the Visibility Window 

% Define stations and their properties
stations = struct( ...
    'name', {'KOUROU', 'TROLL', 'SVALBARD'}, ...
    'min_elevation', [deg2rad(6), deg2rad(0), deg2rad(8)], ...
    'time_grid', {'60s', '30s', '60s'} ...
);

% Loop through each station for computations
for s = 1:3 % Iterate over all stations
    station = stations(s);
    station_name = station.name;
    min_elevation = station.min_elevation(s);  % Access the min_elevation directly
    time_grid = station.time_grid;

    sigma_elevation = deg2rad(125*1e-3);  %[deg]
    sigma_azimuth = deg2rad(125*1e-3);  %[deg]
    sigma_range = 0.01;       %[km]
    cov = diag([sigma_range^2 sigma_azimuth^2 sigma_elevation^2]);

    % Select appropriate time span and states based on the grid
    if strcmp(time_grid, '60s')
        tspan = tspan_visibility_60s;
        reci = reci_60s;
        veci = veci_60s;
        xx = xx_60s;
    elseif strcmp(time_grid, '30s')
        tspan = tspan_visibility_30s;
        reci = reci_30s;
        veci = veci_30s;
        xx = xx_30s;
    end

    % Compute SGP4 measurements
    range_SGP4 = zeros(size(tspan));
    azimuth_SGP4 = zeros(size(tspan));
    elevation_SGP4 = zeros(size(tspan));

     % Compute Keplerian measurements
    range_kep = zeros(size(tspan));
    azimuth_kep = zeros(size(tspan));
    elevation_kep = zeros(size(tspan));

    % Initilize SPG4 measurements noise
    range_SGP4_noise = zeros(size(tspan));
    elevation_SGP4_noise= zeros(size(tspan));
    azimuth_SGP4_noise = zeros(size(tspan));

    for i = 1:length(tspan)
        [range_SGP4(i), azimuth_SGP4(i), elevation_SGP4(i)] = get_rae(station_name, tspan(i), [reci(i, :) veci(i, :)]);
        [range_kep(i), azimuth_kep(i), elevation_kep(i)] = get_rae(station_name, tspan(i), xx(i, :));

        % Add noise to measurements
        temp = mvnrnd([range_SGP4(i); elevation_SGP4(i); azimuth_SGP4(i)], cov); 
        range_SGP4_noise(i) = temp(1);  % Noisy range
        elevation_SGP4_noise(i) = temp(2);  % Noisy elevation
        azimuth_SGP4_noise(i) = temp(3);  % Noisy azimuth

    end
    
    % Filter based on minimum elevation
    valid_indices_SGP4 = elevation_SGP4 >= min_elevation;
    valid_indices_kep = elevation_kep >= min_elevation;
    valid_indices_SGP4_noise = elevation_SGP4_noise >= min_elevation;

    if s == 1
       three.visindex_kourou = valid_indices_SGP4_noise;
       three.measurements_KOUROU = [range_SGP4_noise;  azimuth_SGP4_noise; elevation_SGP4_noise];
    elseif s == 2
        three.visindex_troll = valid_indices_SGP4_noise;
       three.measurements_TROLL = [range_SGP4_noise;  azimuth_SGP4_noise; elevation_SGP4_noise];
    elseif s == 3
        three.visindex_svalbard = valid_indices_SGP4_noise;
        three.measurements_SVALBARD = [range_SGP4_noise;  azimuth_SGP4_noise; elevation_SGP4_noise];
    end
    
    % Filter based on minimum elevation

    % three.visibility_indexes{s} = valid_indices_SGP4;
    % 
    % three.range{s} = range_SGP4(three.visibility_indexes{s});
    % three.azimuth{s} = azimuth_SGP4(three.visibility_indexes{s});
    % three.elevation{s} = elevation_SGP4(three.visibility_indexes{s});

    % Convert ET to UTC and datenum for plotting
    utc_dates_SGP4 = arrayfun(@(t) cspice_et2utc(t, 'ISOC', 3), tspan, 'UniformOutput', false);
    utc_dates_kep = arrayfun(@(t) cspice_et2utc(t, 'ISOC', 3), tspan, 'UniformOutput', false);

    date_nums_SGP4 = datenum(utc_dates_SGP4, 'yyyy-mm-ddTHH:MM:SS.FFF');
    date_nums_kep = datenum(utc_dates_kep, 'yyyy-mm-ddTHH:MM:SS.FFF');

    if any(valid_indices_kep)
    % Get the first and last valid times
    start_time = date_nums_kep(find(valid_indices_kep, 1, 'first')); % First visible time
    end_time = date_nums_kep(find(valid_indices_kep, 1, 'last'));   % Last visible time
    
    % Display the visibility period
    fprintf('%s visibility windows [UTC]: %s to %s \n\n', station_name, datestr(start_time, 'dd-mm-yyyy HH:MM:SS'), datestr(end_time, 'dd-mm-yyyy HH:MM:SS'));
    else
        % No visibility found
        fprintf('%s visibility windows [UTC]: None during the specified interval.\n\n', station_name);
    end

    % Create the figure for elevation plots
    figure
    scatter(rad2deg(azimuth_kep(valid_indices_kep)), rad2deg(elevation_kep(valid_indices_kep)), 150, 'red', 'Marker', 'x');
    hold on;
    scatter(rad2deg(azimuth_SGP4(valid_indices_SGP4)), rad2deg(elevation_SGP4(valid_indices_SGP4)), 150, 'blue', 'Marker', 'x');
    scatter(rad2deg(azimuth_SGP4_noise(valid_indices_SGP4_noise)), rad2deg(elevation_SGP4_noise(valid_indices_SGP4_noise)), 150, 'green', 'Marker', 'x');
    title(['Azimuth - Elevation @', station_name]);
    xlabel('Azimuth [deg]');
    ylabel('Elevation [deg]');
    ylim([rad2deg(min_elevation), rad2deg(max(elevation_SGP4)) + 2]);
    legend('Keplerian', 'SGP4', 'Noise', 'Location', 'Best');
    grid on;

end

%% 3) Solve The Navigation Problem 
close all; clc;

initial_state = [reci1(end, :) veci1(end, :)];

sigma_elevation = deg2rad(125*1e-3);  %[deg]
sigma_azimuth = deg2rad(125*1e-3);  %[deg]
sigma_range = 0.01;  %[km]
three.cost = [30 35 35];

sigma_meas = diag([sigma_range^2; sigma_azimuth^2; sigma_elevation^2]);
%W_m = diag(1./sigma_meas);
W_m = inv(sqrtm(sigma_meas));

three.stations = {'KOUROU', 'TROLL', 'SVALBARD'};
three.min_elevation = [deg2rad(6), deg2rad(0), deg2rad(8)];

% Point a - Measurements from only Kourou and Keplerian Motion

options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display', 'off');

[x_a,resnorm_a,residual_a,~,~,~,jacobian_a] = lsqnonlin(@(x) costfunction(x, W_m, three, 'a'), initial_state, [], [], options);
   
position_error_a = norm(reci1(end, 1:3) - x_a(1:3));
velocity_error_a = norm(veci1(end, 1:3) - x_a(4:6));

jacobian_a = full(jacobian_a);
P_a = resnorm_a/(length(residual_a)-length(x_a)).*inv(jacobian_a.'*jacobian_a);
postrace_a = sqrt(trace(P_a(1:3,1:3)));
veltrace_a = sqrt(trace(P_a(4:6,4:6)));

T = PcarToPkep(x_a, mu);

P_a_kep = T*P_a*T';

% Extract the standard deviations (square root of the diagonal elements)
sigma_a = sqrt(P_a_kep(1, 1));  % Standard deviation of semimajor axis
sigma_i = sqrt(P_a_kep(2, 2));  % Standard deviation of inclination

% Display the results
fprintf('Kourou Measurements & Keplerian Propagation\n');
fprintf('-------------------------------------------\n');
% Estimated State
fprintf('Estimated Position [km]:     [%8.3f %8.3f %8.3f]\n', x_a(1:3));
fprintf('Estimated Velocity [km/s]:   [%8.5f %8.5f %8.5f]\n', x_a(4:6));
fprintf('Position Error Norm: %.3f km\n', position_error_a);
fprintf('Velocity Error Norm: %.5f km/s\n', velocity_error_a);
fprintf('Position Covariance Trace: %.5f km\n', postrace_a);
fprintf('Velocity Covariance Trace: %.5f km/s\n', veltrace_a);
fprintf('Standard deviation of semimajor axis:  %.5f km\n', sigma_a);
fprintf('Standard deviation of inclination:  %.6f deg\n', rad2deg(sigma_i));
fprintf('-------------------------------------------\n\n');

% Point b -  Measurements from every station and Keplerian Motion

options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display', 'off');

[x_b,resnorm_b,residual_b,~,~,~,jacobian_b] = lsqnonlin(@(x) costfunction(x, W_m, three, 'b'), initial_state, [], [], options);
   
position_error_b = norm(reci1(end, 1:3) - x_b(1:3));
velocity_error_b = norm(veci1(end, 1:3) - x_b(4:6));

jacobian_b = full(jacobian_b);
P_b = resnorm_b/(length(residual_b)-length(x_b)).*inv(jacobian_b.'*jacobian_b);

postrace_b = sqrt(trace(P_b(1:3,1:3)));
veltrace_b = sqrt(trace(P_b(4:6,4:6)));

T = PcarToPkep(x_b, mu);

P_b_kep = T*P_b*T';

% Extract the standard deviations (square root of the diagonal elements)
sigma_a = sqrt(P_b_kep(1, 1));  % Standard deviation of semimajor axis
sigma_i = sqrt(P_b_kep(2, 2));  % Standard deviation of inclination

% Display the results
fprintf('All Stations Measurements & Keplerian Propagation\n');
fprintf('-------------------------------------------------\n');
% Estimated State
fprintf('Estimated Position [km]:     [%8.3f %8.3f %8.3f]\n', x_b(1:3));
fprintf('Estimated Velocity [km/s]:   [%8.5f %8.5f %8.5f]\n', x_b(4:6));
fprintf('Position Error Norm: %.3f km\n', position_error_b);
fprintf('Velocity Error Norm: %.5f km/s\n', velocity_error_b);
fprintf('Position Covariance Trace: %.5f km\n', postrace_b);
fprintf('Velocity Covariance Trace: %.5f km/s\n', veltrace_b);
fprintf('Standard deviation of semimajor axis:  %.5f km\n', sigma_a);
fprintf('Standard deviation of inclination:  %.6f deg\n', rad2deg(sigma_i));
fprintf('-------------------------------------------------\n\n');

% Point c -  Measurements from every station and J2 Motion

options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display', 'off');

[x_c,resnorm_c,residual_c,exitflag,~,~,jacobian_c] = lsqnonlin(@(x) costfunction(x, W_m, three, 'c'), initial_state, [], [], options);
   
position_error_c = norm(reci1(end, 1:3) - x_c(1:3));
velocity_error_c = norm(veci1(end, 1:3) - x_c(4:6));

jacobian_c = full(jacobian_c);
P_c = resnorm_c/(length(residual_c)-length(x_c)).*inv(jacobian_c.'*jacobian_c);

postrace_c = sqrt(trace(P_c(1:3,1:3)));
veltrace_c = sqrt(trace(P_c(4:6,4:6)));

T = PcarToPkep(x_c, mu);

P_c_kep = T*P_c*T';

% Extract the standard deviations (square root of the diagonal elements)
sigma_a = sqrt(P_c_kep(1, 1));  % Standard deviation of semimajor axis
sigma_i = sqrt(P_c_kep(2, 2));  % Standard deviation of inclination

% Display the results
fprintf('All Stations Measurements & J2 Propagation\n');
fprintf('------------------------------------------\n');
fprintf('Estimated Position [km]:     [%8.3f %8.3f %8.3f]\n', x_c(1:3));
fprintf('Estimated Velocity [km/s]:   [%8.5f %8.5f %8.5f]\n', x_c(4:6));
fprintf('Position Error Norm: %.3f km\n', position_error_c);
fprintf('Velocity Error Norm: %.5f km/s\n', velocity_error_c);
fprintf('Position Covariance Trace: %.4f km\n', postrace_c);
fprintf('Velocity Covariance Trace: %.7f km/s\n', veltrace_c);
fprintf('Standard deviation of semimajor axis:  %.4f km\n', sigma_a);
fprintf('Standard deviation of inclination:  %.7f deg\n', rad2deg(sigma_i));
fprintf('------------------------------------------\n\n');


%% 4) Trade Off Analysis
clc;
options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display', 'off');

station_combinations = nchoosek(1:3, 2);
three.cost = [30 35 35];

% Initialize variables to store results
best_combination = [];
min_score = Inf;  % Start with a very high score
results = [];

% Weights for the combined objective (adjust as needed)
w_c = 1;  % Weight for cost
w_a = 100;  % Weight for semi-major axis 
w_i = 10^5; % Weight for inclination
 
for i = 1:length(station_combinations)
       three.combinations = station_combinations(i, :);
       [x_d,resnorm_d,residual_d,exitflag,~,~,jacobian_d] = lsqnonlin(@(x) costfunction(x, W_m, three, 'd'), initial_state, [], [], options);
           
        % position_error_d = norm(reci1(end, 1:3) - x_d(1:3));
        % velocity_error_d = norm(veci1(end, 1:3) - x_d(4:6));
        % fprintf('Point c - Position Error Norm: %.5f km\n', position_error_d);
        % fprintf('Point c - Velocity Error Norm: %.5f km/s\n\n', velocity_error_d);
        
        jacobian_d = full(jacobian_d);
        P_d = resnorm_d/(length(residual_d)-length(x_d)).*inv(jacobian_d.'*jacobian_d);
        
        postrace_d = sqrt(trace(P_d(1:3,1:3)));
        veltrace_d = sqrt(trace(P_d(4:6,4:6)));
        
        T = PcarToPkep(x_d, mu);
        
        P_d_kep = T*P_d*T';
        
        % Extract the standard deviations (square root of the diagonal elements)
        sigma_a = sqrt(P_d_kep(1, 1));  % Standard deviation of semimajor axis
        sigma_i = sqrt(P_d_kep(2, 2));  % Standard deviation of inclination
        
        total_cost = sum(three.cost(three.combinations));
   
        % Compute the combined score
        score = w_c * total_cost + w_a * sigma_a + w_i * sigma_i;
        
        % Display the results
        fprintf('Selected Stations: %s\n', strjoin(three.stations(three.combinations), ', '));
        fprintf('Standard deviation of semimajor axis (sigma_a):  %.5f km\n', sigma_a);
        fprintf('Standard deviation of inclination (sigma_i):  %.6f deg\n', rad2deg(sigma_i));
        fprintf('Total Cost: %.0fK\n', total_cost);
        fprintf('SCORE: %.0f\n\n', score) % Winner is the one with lowest score (like in golf)             
   end

%% 5) Long-Term Visibility Analysis 
close all;

% Define 72h window from start of visibility
et_long_start = et_0;
et_long_end = et_0 + 3*86400;  % 72 hours in ET seconds

% Native time grids per station
three_long.tspan{1} = et_long_start:60:et_long_end;   % KOUROU (60s)
three_long.tspan{2} = et_long_start:30:et_long_end;   % TROLL (30s)
three_long.tspan{3} = et_long_start:60:et_long_end;   % SVALBARD (60s)

% Storage
reci_sgp4 = cell(1, 3);
veci_sgp4 = cell(1, 3);

% SGP4 Propagation
for s = 1:3
    tspan = three_long.tspan{s};
    reci = zeros(length(tspan), 3);
    veci = zeros(length(tspan), 3);
    
    for i = 1:length(tspan)
        tsince = (tspan(i) - et_ref)/60.0;  % minutes since TLE epoch
        [~, rteme, vteme] = sgp4(satrec, tsince);
        ttt = cspice_unitim(tspan(i), 'ET', 'TDT') / cspice_jyear() / 100;
        [reci_i, veci_i] = teme2eci(rteme, vteme, [0;0;0], ttt, ddpsi, ddeps);
        reci(i,:) = reci_i';
        veci(i,:) = veci_i';
    end

    reci_sgp4{s} = reci;
    veci_sgp4{s} = veci;
end

% Station setup
stations = {'KOUROU', 'TROLL', 'SVALBARD'};
min_el = [deg2rad(6), deg2rad(0), deg2rad(8)];
az_deg = cell(1, 3);
el_deg = cell(1, 3);
tplot = cell(1, 3);

for s = 1:3
    tspan = three_long.tspan{s};
    reci = reci_sgp4{s};
    veci = veci_sgp4{s};
    
    az = zeros(size(tspan));
    el = zeros(size(tspan));
    
    for i = 1:length(tspan)
        [~, az(i), el(i)] = get_rae(stations{s}, tspan(i), [reci(i,:) veci(i,:)]);
    end

    % --- PASS ANALYSIS (using full time grid, not just visible samples) ---
    pass_flags = el > min_el(s);  % Logical: 1 if elevation > threshold
    pass_edges = diff([0; pass_flags(:); 0]);
    pass_starts = find(pass_edges == 1);
    pass_ends = find(pass_edges == -1) - 1;
    num_passes = length(pass_starts);

    durations = zeros(1, num_passes);
    max_elev = zeros(1, num_passes);

    for p = 1:num_passes
        idx = pass_starts(p):pass_ends(p);
        durations(p) = (tspan(idx(end)) - tspan(idx(1)));  % seconds
    
        % Convert to readable date format
        start_dn = datenum(cspice_et2utc(tspan(idx(1)), 'ISOC', 3), 'yyyy-mm-ddTHH:MM:SS.FFF');
        end_dn = datenum(cspice_et2utc(tspan(idx(end)), 'ISOC', 3), 'yyyy-mm-ddTHH:MM:SS.FFF');
    
        start_str = datestr(start_dn, 'dd-mmm-yyyy HH:MM:SS');
        end_str = datestr(end_dn, 'dd-mmm-yyyy HH:MM:SS');
    
        fprintf('Pass %02d: %s --> %s | Duration: %.1f min\n', p, start_str, end_str, durations(p)/60);
    end

    % Convert for scatter plot
    visible = el > min_el(s);
    az_deg{s} = rad2deg(az(visible));
    el_deg{s} = rad2deg(el(visible));
    utc_strs = arrayfun(@(et) cspice_et2utc(et, 'ISOC', 3), tspan(visible), 'UniformOutput', false);
    tplot{s} = datetime(utc_strs, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');

    % Output
    fprintf('--- %s VISIBILITY SUMMARY (72h) ---\n', stations{s});
    fprintf('Number of Passes: %d\n', num_passes);
    fprintf('Average Pass Duration: %.1f min\n', mean(durations)/60);
    
    % Optional: store if needed
    visibility_stats(s).station = stations{s};
    visibility_stats(s).num_passes = num_passes;
    visibility_stats(s).durations = durations;
end

%% Plot Elevation vs Time for each station
for s = 1:3
    figure
    scatter(tplot{s}, el_deg{s}, 'red', 'Marker', '*');
    title(['Elevation over 72h (SGP4) @ ', stations{s}]);
    xlabel('Date');
    ylabel('Elevation [deg]');
    ylim([rad2deg(min_el(s)), 90]);
    grid on;
    datetick('x', 'dd-mm-yyyy HH:MM', 'keepticks');
end


% Save Figures

% saveAllFiguresToPath('/Users/pierazhi/Documents/MATLAB/GNC/Assignment 2 - Ex2/Images')

%% Functions

function [range, azimuth, elevation] = get_rae(station, et, xx)
%--------------------------------------------------------------------------
% Function: get_rae
% Author: Piercarlo Fontana
% Date: 30 May 2025
%
% Computes the range, azimuth, and elevation of a satellite with respect 
% to a ground station at a given epoch time, based on the satellite's 
% Earth-centered inertial (ECI) state.
%
% Inputs:
%   station   - String name of the ground station (e.g., 'KOUROU')
%   et        - Ephemeris Time (ET), in seconds past J2000 TDB
%   xx        - Satellite state vector [x, y, z, vx, vy, vz] in km and km/s
%
% Outputs:
%   range     - Slant range from station to satellite [km]
%   azimuth   - Azimuth angle relative to local ENU frame [rad]
%   elevation - Elevation angle above local horizon [rad]
%--------------------------------------------------------------------------

% Compute station position in ECI
rv_station_eci = cspice_spkezr(station, et, 'J2000', 'NONE', 'Earth');

% Compute station-satellite vector in ECI
rv_station_sat_eci = [xx(1, 1:3) xx(1, 4:6)]' - rv_station_eci;

% Define station name
topoFrame = [station, '_TOPO'];

% Transformation from ECI to topocentric frame
ROT_ECI2TOPO = cspice_sxform('J2000', topoFrame, et);

% Convert station-satellite into topocentric frame
rv_station_sat_topo = ROT_ECI2TOPO * rv_station_sat_eci;

% Compute range, azimuth and elevation using cspice_xfmsta
rll_station_sat = cspice_xfmsta(rv_station_sat_topo, 'RECTANGULAR', 'LATITUDINAL', 'Earth');

range     = rll_station_sat(1);   % [km]
azimuth   = rll_station_sat(2);   % [rad]
elevation = rll_station_sat(3);   % [rad]

end

function T = PcarToPkep(state, mu_value)
%--------------------------------------------------------------------------
% Function: PcarToPkep
% Author: Piercarlo Fontana
% Date: 30 May 2025
%
% Computes the Jacobian matrix of the transformation from Cartesian 
% coordinates to selected Keplerian orbital elements: semimajor axis and 
% inclination.
%
% The Jacobian is derived symbolically and evaluated at the specified state.
%
% Inputs:
%   state     - Cartesian state vector [x; y; z; vx; vy; vz] (6x1)
%   mu_value  - Gravitational parameter [km^3/s^2]
%
% Output:
%   T         - Jacobian matrix (2x6) of [a; i] with respect to [x; y; z; vx; vy; vz]
%--------------------------------------------------------------------------

    % Define symbolic variables for position, velocity, and gravitational parameter
    syms x y z vx vy vz real
    mu = sym('mu', 'real');

    % Define position and velocity vectors
    R = [x; y; z];  % Position vector
    V = [vx; vy; vz];  % Velocity vector
    
    % Compute magnitudes (norms)
    r = sqrt(x^2 + y^2 + z^2);  % Position magnitude (distance from the origin)
    v = sqrt(vx^2 + vy^2 + vz^2);  % Velocity magnitude
    
    % Specific orbital energy (epsilon)
    epsilon = (v^2) / 2 - mu / r;
    
    % Semimajor axis (a) based on orbital energy
    a = -mu / (2 * epsilon);  % Semimajor axis
    
    % Inclination (i)
    H = cross(R, V);  % Angular momentum vector
    h = norm(H);  % Magnitude of angular momentum
    i = acos(H(3) / h);  % Inclination

    % Compute the Jacobian matrix for semimajor axis (a) and inclination (i)
    T_sym = jacobian([a; i], [x, y, z, vx, vy, vz]);

    % Substitute numerical values into the Jacobian matrix
    T = subs(T_sym, {x, y, z, vx, vy, vz, mu}, {state(1), state(2), state(3), state(4), state(5), state(6), mu_value});

    % Convert the symbolic Jacobian to a numerical matrix
    T = double(T);  % Evaluate the symbolic Jacobian matrix numerically

end

function [residual, index, index_pert] = costfunction(x, W_m, three, caso)
%--------------------------------------------------------------------------
% Function: costfunction
% Author: Piercarlo Fontana
% Date: 30 May 2025
%
% Computes the residuals between predicted and observed measurements 
% (range, azimuth, elevation) for various tracking station configurations.
%
% Supports multiple scenarios based on the input case:
%   'a' - Kourou only, unperturbed Keplerian propagation
%   'b' - All stations, unperturbed Keplerian propagation
%   'c' - All stations, J2-perturbed propagation
%   'd' - Station subset combination (e.g., trade-off cases)
%
% Inputs:
%   x      - Estimated state vector
%   W_m    - Weighting matrix for residual computation
%   three  - Struct containing station data, measurements, and visibility indices
%   caso   - Case type identifier: 'a', 'b', 'c', or 'd'
%
% Outputs:
%   residual     - Vector of weighted residuals
%   index        - Active indices for visible measurements (case-dependent)
%   index_pert   - Perturbed visibility indices (not used in all cases)
%--------------------------------------------------------------------------
    % Initialize residual as empty
    residual = [];

    % Check the case
if isequal(caso, 'a')
        % Select the station and associated parameters
        s = 1;  % Assuming 'KOUROU'
        station = three.stations{s};
        tspan = three.tspan{s};
        min_elevation = three.min_elevation(s);

        % Initialize predicted measurements
        elevation_pred = zeros(1, length(tspan));
        azimuth_pred = zeros(1, length(tspan));
        range_pred = zeros(1, length(tspan));

        % Propagate state vector to measurement epochs using Keplerian motion
        [~, ~, ~, xx] = keplerian_propagator(tspan, x, 'Earth', 'keplerian');

        % Compute predicted measurements (RAE) at each epoch
        for i = 1:length(tspan)
            [range_pred(i), azimuth_pred(i), elevation_pred(i)] = get_rae(station, tspan(i), xx(i, :));
        end

        % Identify valid predictions based on elevation cutoff
        index = three.visindex_kourou;

        % Extract valid predicted measurements
        meas_real = three.measurements_KOUROU(:, index);

        % Compute residuals
        for k = 1:length(meas_real)
            % Weighted residual using angular difference
            %diff_meas_weighted = W_m * angdiff(meas_pred(:, k), meas_real(:, k));

            diff_azimuth = angdiff(azimuth_pred(index), meas_real(2, :));
            diff_elevation = angdiff(elevation_pred(index), meas_real(3, :));
            diff_range = range_pred(index) - meas_real(1, :);

            %residual = [residual; diff_meas_weighted(:)];  % Append to residuals
            diff_overall = W_m*[diff_range; diff_elevation; diff_azimuth];
            residual = diff_overall';
        end

elseif isequal(caso, 'b')
            for s = 1:3
                % Select the station and associated parameters
                station = three.stations{s};
                tspan = three.tspan{s};
                min_elevation = three.min_elevation(s);
                
                % Initialize predicted measurements
                elevation_pred = zeros(1, length(tspan));
                azimuth_pred = zeros(1, length(tspan));
                range_pred = zeros(1, length(tspan));
                
                % Propagate state vector to measurement epochs using Keplerian motion
                [~, ~, ~, xx] = keplerian_propagator(tspan, x, 'Earth', 'keplerian');
                
                % Compute predicted measurements (RAE) at each epoch
                for i = 1:length(tspan)
                    [range_pred(i), azimuth_pred(i), elevation_pred(i)] = get_rae(station, tspan(i), xx(i, :));
                end
                
                if s == 1
                   index = three.visindex_kourou;
                   meas_real = three.measurements_KOUROU(:, index);
                elseif s == 2
                   index = three.visindex_troll;
                   meas_real = three.measurements_TROLL(:, index);                
                elseif s == 3
                   index = three.visindex_svalbard; 
                   meas_real = three.measurements_SVALBARD(:, index);                
                end
                
                % Compute residuals (weighted by W_m) for this station
                for k = 1:length(meas_real)
                    % Compute the residual for the current station
                    
                    diff_range = range_pred(index) - meas_real(1, :);
                    diff_azimuth = angdiff(azimuth_pred(index), meas_real(2, :));
                    diff_elevation = angdiff(elevation_pred(index), meas_real(3, :));

                    diff_overall = W_m*[diff_range; diff_elevation; diff_azimuth];
                    residual = [residual; diff_overall(:)];
                end
            end
 elseif isequal(caso, 'c')
            for s = 1:3
                % Select the station and associated parameters
                station = three.stations{s};
                tspan = three.tspan{s};
                min_elevation = three.min_elevation(s);
                
                % Initialize predicted measurements
                elevation_pred = zeros(1, length(tspan));
                azimuth_pred = zeros(1, length(tspan));
                range_pred = zeros(1, length(tspan));
                
                % Propagate state vector to measurement epochs using Keplerian motion
                [~, ~, ~, xx] = keplerian_propagator(tspan, x, 'Earth', 'J2');
                
                % Compute predicted measurements (RAE) at each epoch
                for i = 1:length(tspan)
                    [range_pred(i), azimuth_pred(i), elevation_pred(i)] = get_rae(station, tspan(i), xx(i, :));
                end
                
                 if s == 1
                   index = three.visindex_kourou;
                   meas_real = three.measurements_KOUROU(:, index);
                elseif s == 2
                   index = three.visindex_troll;
                   meas_real = three.measurements_TROLL(:, index);                
                elseif s == 3
                   index = three.visindex_svalbard; 
                   meas_real = three.measurements_SVALBARD(:, index);                
                 end
                
                % Compute residuals (weighted by W_m) for this station
                for k = 1:length(meas_real)
                    % Compute the residual for the current station
                    
                    diff_range = range_pred(index) - meas_real(1, :);
                    diff_azimuth = angdiff(azimuth_pred(index), meas_real(2, :));
                    diff_elevation = angdiff(elevation_pred(index), meas_real(3, :));

                    diff_overall = W_m*[diff_range; diff_elevation; diff_azimuth];
                    residual = [residual; diff_overall(:)];

                end
            end
        elseif isequal(caso, 'd')
            for s = three.combinations
                % Select the station and associated parameters
                station = three.stations{s};
                tspan = three.tspan{s};
                min_elevation = three.min_elevation(s);
                
                % Initialize predicted measurements
                elevation_pred = zeros(1, length(tspan));
                azimuth_pred = zeros(1, length(tspan));
                range_pred = zeros(1, length(tspan));
                
                % Propagate state vector to measurement epochs using Keplerian motion
                [~, ~, ~, xx] = keplerian_propagator(tspan, x, 'Earth', 'J2');
                
                % Compute predicted measurements (RAE) at each epoch
                for i = 1:length(tspan)
                    [range_pred(i), azimuth_pred(i), elevation_pred(i)] = get_rae(station, tspan(i), xx(i, :));
                end
                
                 if s == 1
                   index = three.visindex_kourou;
                   meas_real = three.measurements_KOUROU(:, index);
                elseif s == 2
                   index = three.visindex_troll;
                   meas_real = three.measurements_TROLL(:, index);                
                elseif s == 3
                   index = three.visindex_svalbard; 
                   meas_real = three.measurements_SVALBARD(:, index);                
                 end
                  
                % Compute residuals (weighted by W_m) for this station
                for k = 1:length(meas_real)
                    % Compute the residual for the current station
                    
                    diff_range = range_pred(index) - meas_real(1, :);
                    diff_azimuth = angdiff(azimuth_pred(index), meas_real(2, :));
                    diff_elevation = angdiff(elevation_pred(index), meas_real(3, :));

                    diff_overall = W_m*[diff_range; diff_elevation; diff_azimuth];
                    residual = [residual; diff_overall(:)];
                end
            end
        end
   residual = residual(:);
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