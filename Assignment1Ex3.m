% Spacecraft Guidance & Navigation
% Assignment # 1, Exercise 3
% Author: Piercarlo Fontana

%% 1) Plot Debris Density
clearvars; clc; close all; set(0,'DefaultFigureVisible','on')
format long G

const.Re = 6378.1366;
const.hi = 800;
const.hf = 1000;
const.mu_d = 398600.435;
const.DU = 7178.1366;
const.TU = sqrt(const.DU^3 / const.mu_d);
const.MU = 1000;
const.VU = const.DU/const.TU;

const.ri = [const.Re + const.hi, 0, 0]'/const.DU;
const.vi = [0, sqrt(const.mu_d / (const.hi + const.Re)), 0]'/const.VU;

const.rf = [const.Re + const.hf, 0, 0]'/const.DU;
const.vf = [0, sqrt(const.mu_d / (const.hf + const.Re)) * cosd(0.75), sqrt(const.mu_d / (const.hf + const.Re)) * sind(0.75)]'/const.VU;

const.Isp = 3120 / const.TU;
const.Tmax = 3 * 10^-3 * const.TU^2 / (const.MU * const.DU);
const.g0 = 9.81 * 10^-3 * const.TU^2 / const.DU;
const.mu = 1;
const.rho0 = (750 + const.Re)/const.DU;
const.k1 = 10^-5;
const.k2 = 10^-4;   

const.initial_state = [const.ri; const.vi; 1];

density = @(rho) const.k1 ./ ( const.k2 + ((rho - 750)/const.DU).^2);
h = linspace(700, 1100, 10000);

figure
plot(h, density(h), 'b-', 'LineWidth', 2)
% Add labels and title for clarity
xlabel('Distance from the Earth Center [km]')
ylabel('Debris Spatial Density')
title('Spatial Density of Space Debris')
grid on

% For the circular orbit (equatorial, no inclination)
T1 = (2*pi*sqrt((const.Re + const.hi)^3/const.mu_d));
tspan1 = linspace(0, T1);
options = odeset( 'RelTol', 1e-12, 'AbsTol', 1e-12);
[~, Y_equa] = ode113(@(t,y) ode_j2drag(t, y, const.mu_d, 'no'), tspan1, [const.Re + const.hi, 0, 0, 0, sqrt(const.mu_d / (const.hi + const.Re)), 0], options);

% For the circular orbit (equatorial, inclination)

T2 = (2*pi*sqrt((const.Re + const.hf)^3/const.mu_d));
tspan2 = linspace(0, T2);
options = odeset( 'RelTol', 1e-12, 'AbsTol', 1e-12);
[~, Y_incl] = ode113(@(t,y) ode_j2drag(t, y, const.mu_d, 'no'), tspan2, [const.Re + const.hf, 0, 0, 0, sqrt(const.mu_d / (const.hf + const.Re)) * cosd(0.75), sqrt(const.mu_d / (const.hf + const.Re)) * sind(0.75)], options);

%Plot the orbits
figure
%plotPlanet(3, [0 0 0], gca, 1)
hold on
plot3(Y_equa(1:10, 1), Y_equa(1:10, 2), Y_equa(1:10, 3), 'b-', 'LineWidth', 2) % Circular orbit in blue
plot3(Y_incl(1:10, 1), Y_incl(1:10, 2), Y_incl(1:10, 3), 'r-', 'LineWidth', 2) % Circular orbit in red
% Mark the initial and final points for both orbits
plot3(Y_equa(1, 1), Y_equa(1, 2), Y_equa(1, 3), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');  % Initial point (green)
plot3(Y_incl(1, 1), Y_incl(1, 2), Y_incl(1, 3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % Final point (blue)
plot3(Y_equa(end-10:end, 1), Y_equa(end-10:end, 2), Y_equa(end-10:end, 3), 'b-', 'LineWidth', 2) % Circular orbit in blue
plot3(Y_incl(end-10:end, 1), Y_incl(end-10:end, 2), Y_incl(end-10:end, 3), 'r-', 'LineWidth', 2) % Circular orbit in red
% Labels and title
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
title('Initial and Final Orbit @ ECI)')
grid on
view(3)
legend('Initial Orbit', 'Final Orbit', 'Initial Point', 'Final Point', 'Location', 'best')

fprintf('\n--- Initial Orbit  ---\n');
fprintf('x0    = %12.6f \n', const.ri(1)*const.DU);
fprintf('y0    = %12.6f \n', const.ri(2)*const.DU);
fprintf('z0    = %12.6f \n', const.ri(3)*const.DU);
fprintf('x0dot = %12.6f \n', const.vi(1)*const.VU);
fprintf('y0dot = %12.6f \n', const.vi(2)*const.VU);
fprintf('z0dot = %12.6f \n', const.vi(3)*const.VU);

fprintf('\n--- Final Orbit  ---\n');
fprintf('x0    = %12.6f \n', const.rf(1)*const.DU);
fprintf('y0    = %12.6f \n', const.rf(2)*const.DU);
fprintf('z0    = %12.6f \n', const.rf(3)*const.DU);
fprintf('x0dot = %12.6f \n', const.vf(1)*const.VU);
fprintf('y0dot = %12.6f \n', const.vf(2)*const.VU);
fprintf('z0dot = %12.6f \n', const.vf(3)*const.VU);

% Do Not Run This Section, as the solution is not always found (it takes time); instead, run the following point for the already found costate solution

close all; clearvars -except const

run_this_section = false;  % Set to true to run the block

if run_this_section
    const.Tmax = 3 * 10^-3 * const.TU^2 / (const.MU * const.DU);

    %Initial Guess
    
    sol_converged = false;
    
    options = optimoptions('fsolve', 'Display', 'iter-detailed', 'MaxFunctionEvaluations', 3000);
    max_iter = 100;
    iter = 1;
    
    while ~sol_converged && iter < max_iter
        % Update initial guess for lambda0 with a random component
        lambda0 = [randi([-250, 250], 1, 6)'; randi([0, 250], 1, 1)']*rand(1);
        X0 = [lambda0; 20*pi + 0.5*pi*rand(1)];
    
        [sol, ~, exitflag] = fsolve(@(X) zerofinder(X, const), X0, options);
    
        if exitflag == 1
            sol_converged = true;
        else
            iter = iter + 1;
        end
    end
end

%% 4) Solve the problem
clc; clearvars -except const

sol = [-214.981130561796
-10.3658462643537
0.885592822924761
-10.3928939074737
-214.610360261114
-112.945247816343
2.59644696327156
64.4801061435570];

fprintf('\n--- Costates for T = 3N  ---\n');
fprintf('lambda0_rx    = %12.4f \n', sol(1));
fprintf('lambda0_ry    = %12.4f \n', sol(2));
fprintf('lambda0_rz    = %12.4f \n', sol(3));
fprintf('lambda0_vx    = %12.4f \n', sol(4));
fprintf('lambda0_vy    = %12.4f \n', sol(5));
fprintf('lambda0_vz    = %12.4f \n', sol(6));
fprintf('lambda0_m     = %12.4f \n\n', sol(7));

[errs, YY, errs_km, tspan] = zerofinder(sol, const);

fprintf('Final Mass: %.4f kg\n', YY(end, 7) * const.MU);

fprintf('Final Time: %.4f minutes\n', sol(8) * const.TU / 60);

fprintf('Position Error: %.4f km\n', norm(errs_km(1:3)));
fprintf('Velocity Error: %.4f m/s\n', norm(errs_km(4:6)));

H = zeros(1, length(YY));
primer = zeros(3, length(YY));
primer_NTW = zeros(3, length(YY));

for i = 1:length(YY)
    [H(i), primer(:, i), primer_NTW(:, i)] = Hamiltonian(YY(i, :), const);
end

max_error_fromavg = max(abs(H - movmean(H, round(0.1*length(YY)))));
max_error_fromstart = max(abs(H - H(1)));

figure
hold on
plot(tspan/(2*pi), H, 'b-', 'LineWidth', 1)  % Raw data in blue dots
plot(tspan/(2*pi), movmean(H, round(0.05*length(YY))), 'k-', 'LineWidth', 2)  % Moving average in black
grid on
xlabel('Number of Revolutions')
ylabel('Hamiltonian Value')
title('Hamiltonian Evolution')
legend('Raw Data', 'Moving Average', 'Location', 'best')
hold off

% NTW Components Evolution

figure
hold on
plot(tspan/(2*pi), primer_NTW(1, :), 'r', 'LineWidth', 1.5)
plot(tspan/(2*pi), primer_NTW(2, :), 'g', 'LineWidth', 1.5)
plot(tspan/(2*pi), primer_NTW(3, :), 'b', 'LineWidth', 1.5)
grid on
xlabel('Number of Revolutions')
ylabel('Component')
title('Primer Vector - Components Evolution')
legend('N', 'T', 'W', 'Location', 'best')
hold off

figure
hold on
grid on
plot3(YY(:, 1)*const.DU, YY(:, 2)*const.DU, YY(:, 3)*const.DU, 'b-', 'LineWidth', 1.5)
plot3(const.ri(1)*const.DU, const.ri(2)*const.DU, const.ri(3)*const.DU, 'r.', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Initial Position')
plot3(const.rf(1)*const.DU, const.rf(2)*const.DU, const.rf(3)*const.DU, 'g.', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Final Position')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title('Trajectory for T = 3.000 N @ J2000')
legend('Trajectory', 'Initial Position', 'Final Position', 'Location', 'best')
view(3) % Adjust the angles (azimuth and elevation) as needed
hold off

%% 5) Decrease the thrust

clc; clearvars -except const

const.Re = 6378.1366;
const.hi = 800;
const.hf = 1000;
const.mu_d = 398600.435;
const.DU = 7178.1366;
const.TU = sqrt(const.DU^3 / const.mu_d);
const.MU = 1000;
const.VU = const.DU/const.TU;

const.ri = [const.Re + const.hi, 0, 0]'/const.DU;
const.vi = [0, sqrt(const.mu_d / (const.hi + const.Re)), 0]'/const.VU;

const.rf = [const.Re + const.hf, 0, 0]'/const.DU;
const.vf = [0, sqrt(const.mu_d / (const.hf + const.Re)) * cosd(0.75), sqrt(const.mu_d / (const.hf + const.Re)) * sind(0.75)]'/const.VU;

const.Isp = 3120 / const.TU;
const.Tmax = 3 * 10^-3 * const.TU^2 / (const.MU * const.DU);
const.g0 = 9.81 * 10^-3 * const.TU^2 / const.DU;
const.mu = 1;
const.rho0 = (750 + const.Re)/const.DU;
const.k1 = 10^-5;
const.k2 = 10^-4;   

const.initial_state = [const.ri; const.vi; 1];

N = 10; % Steps to decrease the thrust

Tmax_span = linspace(3 * 10^-3 * const.TU^2 / (const.MU * const.DU),2.860 * 10^-3 * const.TU^2 / (const.MU * const.DU), N);
Tmax_SPAN = Tmax_span/(10^-3*const.TU^2 / (const.MU * const.DU));

sol(:, 1) = [-214.981130561796
            -10.3658462643537
            0.885592822924761
            -10.3928939074737
            -214.610360261114
            -112.945247816343
            2.59644696327156
            64.4801061435570];

options = optimoptions('fsolve', 'Display', 'iter-detailed', 'MaxFunctionEvaluations', 3000);

for i = 2:length(Tmax_span)

    sol_converged = false;

    const.Tmax = Tmax_span(i);
    
    % new initial guess is the solution of the previous iteration
    X = sol(: , i-1);

    % find correct costate initial conditions and propagate the s/c orbit
    while ~sol_converged

            [sol(:, i), ~, exitflag] = fsolve(@(X) zerofinder(X, const), X, options);
        
            if exitflag == 1
                sol_converged = true;
            else
            end
    end
end

ERRORI_KM = zeros(6, size(sol, 2));

for i = 1:size(sol, 2)

    const.Tmax = Tmax_span(i);
    [errs, YY, errs_km(:,i)] = zerofinder(sol(:, i), const);

    fprintf('Tmax = %.5f N\n', Tmax_SPAN(i));
    fprintf('Final Mass: %.4f kg\n', YY(end, 7) * const.MU);
    fprintf('Final Time: %.4f minutes\n', sol(8, i) * const.TU / 60);
    fprintf('Position Error: %.4f km\n', norm(errs_km(1:3, i)));
    fprintf('Velocity Error: %.4f m/s\n', norm(errs_km(4:6, i)));
    fprintf('\n');

    figure
    hold on
    grid on
    plot3(YY(:, 1)*const.DU, YY(:, 2)*const.DU, YY(:, 3)*const.DU, 'b-', 'LineWidth', 1.5)
    plot3(const.ri(1)*const.DU, const.ri(2)*const.DU, const.ri(3)*const.DU, 'r.', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Initial Position')
    plot3(const.rf(1)*const.DU, const.rf(2)*const.DU, const.rf(3)*const.DU, 'g.', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Final Position')
    xlabel('X [DU]')
    ylabel('Y [DU]')
    zlabel('Z [DU]')
    title(sprintf('Trajectory for T = %.3f N @ J2000', Tmax_SPAN(i)))
    legend('Trajectory', 'Initial Position', 'Final Position', 'Location', 'best')
    view(3) % Adjust the angles (azimuth and elevation) as needed
    hold off
end

fprintf('\n--- Costates for T = 2.860 N  ---\n');
fprintf('lambda0_rx    = %12.4f \n', sol(1, end));
fprintf('lambda0_ry    = %12.4f \n', sol(2, end));
fprintf('lambda0_rz    = %12.4f \n', sol(3, end));
fprintf('lambda0_vx    = %12.4f \n', sol(4, end));
fprintf('lambda0_vy    = %12.4f \n', sol(5, end));
fprintf('lambda0_vz    = %12.4f \n', sol(6, end));
fprintf('lambda0_m     = %12.4f \n\n', sol(7, end));

% Save Figures

%saveAllFiguresToPath('/Users/pierazhi/Documents/MATLAB/GNC/Assignment 1 - Ex3/Images')

%% Functions

function [H, primer, primer_NTW] = Hamiltonian(X, const)
%--------------------------------------------------------------------------
% Function: Hamiltonian
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Computes the Hamiltonian value and the primer vector in both Cartesian 
% and NTW (Normal-Tangential-Weak) frames for a spacecraft in optimal 
% control trajectory governed by the Pontryagin Maximum Principle.
%
% Inputs:
%   X     - Augmented state vector [r(3); v(3); m; λ_r(3); λ_v(3); λ_m] (14x1)
%   const - Structure containing constants (Isp, Tmax, g0, mu, rho0, k1, k2)
%
% Outputs:
%   H           - Scalar value of the Hamiltonian
%   primer      - Primer vector (3x1), normalized negative costate of velocity
%   primer_NTW  - Primer vector expressed in NTW frame (3x1)
%--------------------------------------------------------------------------

    Isp = const.Isp;
    Tmax = const.Tmax;
    g0 = const.g0;
    mu = const.mu;
    rho0 = const.rho0;
    k1 = const.k1;
    k2 = const.k2;

    r = X(1:3)';
    v = X(4:6)';
    m = X(7);
    lr = X(8:10)';
    lv = X(11:13)';
    rr = norm(r);
    llv = norm(lv);
    lm = X(14);
    l = [lr;lv;lm];

    q = k1 / (k2 +   (rr - rho0)^2);

    primer = - lv/llv;

    f = [v;
            -(mu/rr^3)*r + (Tmax/m) * primer;
             -Tmax/(Isp*g0)];

    T = v/norm(v);
    W = r/norm(r);
    N = cross(T, W);

    % Transformation matrix from current frame to NTW
    R_NTW = [N'; T'; W'];

    % Transform primer vector to NTW coordinates
    primer_NTW = R_NTW * primer;

    H = q + dot(l, f);

end

function [dy] = ode_j2drag( ~, y, mu, perturbed)
%--------------------------------------------------------------------------
% Function: ode_j2drag
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Computes the time derivative of the state for an orbit under two-body 
% motion, optionally including J2 and atmospheric drag perturbations.
%
% Supports three modes:
%   'no'    - Pure two-body motion (unperturbed)
%   'cart'  - Perturbations modeled in Cartesian coordinates
%   'gauss' - Perturbations modeled using Gauss' variational equations
%
% Inputs:
%   ~         - Time input (ignored)
%   y         - State vector (6x1), either Cartesian or Keplerian elements
%   mu        - Gravitational parameter
%   perturbed - String flag: 'no', 'cart', or 'gauss'
%
% Output:
%   dy - Time derivative of the state (6x1), format depending on mode
%--------------------------------------------------------------------------

% Constants

R_e = astroConstants(23);           
J2 = astroConstants(9);
day = 23*3600+56*60+4;              
we = deg2rad(2*pi/day);            
w_e = [0 0 we]';

if isequal(perturbed, 'no')
    r=y(1:3);
    v=y(4:6);
    rr=norm(r);  
    dy=[v; (-mu/rr^3)*r];

elseif isequal(perturbed, 'cart')
        r=y(1:3);
        v=y(4:6);
        rr=norm(r);

        a_J2_x = (1.5*J2*mu*R_e^2/rr^4) * y(1)/rr*(5*y(3)^2/rr^2-1);
        a_J2_y = (1.5*J2*mu*R_e^2/rr^4) * y(2)/rr*(5*y(3)^2/rr^2-1);
        a_J2_z = (1.5*J2*mu*R_e^2/rr^4) * y(3)/rr*(5*y(3)^2/rr^2-3);

        cd  = 2.1;                                             % Drag coefficient                    [-]
        A2m = 0.0043;                                          % A/m ratio                           [m^2/Kg]
        alt = 1000*(rr-R_e);                                   % Altitude                            [km] 
        rho = atmos(alt);                                      % Air density from Standard Model     [kg/m^3]
        vrel = v - cross(w_e, r);
        Vrel = norm(vrel);                            
        uv = vrel/Vrel;                                        % Relative velocity unit vector
        a_drag_xyz = (-0.5*cd*A2m*rho*(1000*Vrel)^2*uv)/1000;  % drag perturbations                  [km/s^2]
        a_drag_x = a_drag_xyz(1);
        a_drag_y = a_drag_xyz(2);
        a_drag_z = a_drag_xyz(3);
        a_tot_x = a_J2_x + a_drag_x;
        a_tot_y = a_J2_y + a_drag_y;
        a_tot_z = a_J2_z + a_drag_z;
        a_tot_xyz = [a_tot_x a_tot_y a_tot_z]';

        dy=[v; (-mu/rr^3)*r + a_tot_xyz];

elseif isequal(perturbed, 'gauss')
         a  = y(1);   
         e  = y(2);   
         i  = y(3);
         OM = y(4);   
         om = y(5);   
         theta = y(6);
         [r, v] = kep2car(a, e, i, OM, om, theta, mu);
         rr=norm(r);
         a_J2_rsw =J2_rsw(a,e,i,OM,om,theta,mu);
         cd  = 2.1;                                             
         A2m = 0.0043;                                          
         alt = 1000*(rr-R_e);                                   
         rho = atmos(alt);                                     
         vrel = v - cross(w_e, r);                                     
         uv = vrel/norm(vrel);                                        % Relative velocity unit vector
         a_drag_xyz = (-0.5*cd*A2m*rho*(1000*norm(vrel))^2*uv)/1000;  % drag perturbations                  [km/s^2]
         rrr = r';
         vv = v';
         er = rrr/norm(rrr) ;                                   % Tangent unit vector
         ew = cross(rrr,vv)/norm(cross(rrr,vv));                % Normal  unit vector
         es = cross(ew,er); 
         ROT_xyz2rsw = [er;es;ew];
         a_drag_rsw = ROT_xyz2rsw * a_drag_xyz;
         a_tot_rsw = a_J2_rsw + a_drag_rsw;
         a_r = a_tot_rsw(1);  
         a_s = a_tot_rsw(2);  
         a_w = a_tot_rsw(3);

          % Semilatus rectum
          p = a * (1 - e^2);

          % Angular Momentum
          h = sqrt(p * mu);

          % Distance
          r = p / (1 + e * cos(theta));

          % Keplerian element derivatives
          a_dot  = 2 * a^2 * (e * sin(theta) * a_r + p * a_s / r) / h;
          e_dot  = (p * sin(theta) * a_r + a_s * ((p + r) * cos(theta) + r * e)) / h;
          i_dot  = r * cos(theta + om) * a_w / h;
          OM_dot = r * sin(theta + om) * a_w / h / sin(i);
          om_dot  = (-p * cos(theta) * a_r + (p + r) * sin(theta) * a_s) / e / h - r * sin(theta + om) * cos(i) * a_w / h / sin(i);
          theta_dot  = h / r^2 + (p * cos(theta) * a_r - (p + r) * sin(theta) * a_s) / e / h;

          % Derivative vector
          dy = [a_dot, e_dot, i_dot, OM_dot, om_dot, theta_dot]';
end

end

function [dYdt] = rhs_PMP(~,Y, const)
%--------------------------------------------------------------------------
% Function: rhs_PMP
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Defines the right-hand side of the differential equations under the 
% Pontryagin Maximum Principle (PMP) for trajectory optimization with 
% control effort, fuel consumption, and costates.
%
% Inputs:
%   ~     - Time input (ignored)
%   Y     - Augmented state vector [r; v; m; λ_r; λ_v; λ_m] (14x1)
%   const - Structure containing system constants (Isp, Tmax, g0, mu, rho0, k1, k2)
%
% Output:
%   dYdt  - Time derivative of the augmented state vector (14x1)
%--------------------------------------------------------------------------    
        Isp = const.Isp;
        Tmax = const.Tmax;
        g0 = const.g0;
        mu = const.mu;
        rho0 = const.rho0;  
        k1 = const.k1;
        k2 = const.k2;
           
    r = Y(1:3);
    v = Y(4:6);
    m = Y(7);
    l_r = Y(8:10);
    l_v = Y(11:13);
    rr = norm(r);
    ll_v = norm(l_v);

    dYdt = zeros(14,1);

    dYdt(1:3) = v;
    dYdt(4:6) = - mu*r/rr^3 - (Tmax * l_v)/(m*ll_v);
    dYdt(7) = -Tmax/(Isp*g0);
    dYdt(8:10) = - 3*mu * dot(r,l_v)*r / rr^5 + mu*l_v/rr^3 - (2 * k1 * r * (rho0 - rr)) / ((k2 + (rho0 - rr)^2)^2 * rr);
    dYdt(11:13) = -l_r;
    dYdt(14) = -ll_v*Tmax/m^2;

end

function [errs, YY, errs_km, tt] = zerofinder(X, const) 
%--------------------------------------------------------------------------
% Function: zerofinder
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Computes the residuals (errors) of the boundary conditions for a 
% trajectory optimization problem governed by the Pontryagin Maximum Principle.
%
% This function is typically used in conjunction with a root-finding method 
% to enforce boundary constraints by iteratively refining initial guesses.
%
% Inputs:
%   X     - Initial guess for the optimization vector
%   const - Structure with physical and numerical constants
%
% Outputs:
%   errs      - Residual vector (non-dimensional)
%   YY        - Full trajectory matrix
%   errs_km   - Residuals in dimensional units (km/s or km)
%   tt        - Time vector corresponding to YY
%--------------------------------------------------------------------------

        Isp = const.Isp;
        Tmax = const.Tmax;
        g0 = const.g0;
        mu = const.mu;
        rho0 = const.rho0;
        k1 = const.k1;
        k2 = const.k2;
        rf = const.rf;
        vf = const.vf;
        
        tf = X(end);
        
        options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
        [tt,YY] = ode78(@(t, y) rhs_PMP(t,y, const), [0 tf], [const.initial_state; X(1:7)], options);
            
        Yf=YY(end,1:14);
        r = Yf(end,1:3)';
        v = Yf(end,4:6)';
        m = Yf(end,7);
        lr = Yf(end, 8:10)';
        lv = Yf(end, 11:13)';
        rr = norm(r);
        llv = norm(lv);
        lm = Yf(end, 14);
        l = [lr;lv;lm];
        
        q = k1 / (k2 +   (rr - rho0)^2);

        f = [v;
                -(mu/rr^3)*r - (Tmax/m) * lv/llv;
                 -Tmax/(Isp*g0)];

        H = q + dot(l, f);

        errs = [r - rf;
                v - vf;
                lm;
                H];

        errs_km = [errs(1:3)*const.DU; errs(4:6)*const.VU*1000];
           
end

function out = astroConstants(in)

% astroConstants.m - Returns astrodynamic-related physical constants.
%
% PROTOTYPE:
%   out = astro_constants(in)
%
% DESCRIPTION:
%   Returns a row vector of constants, in which there is the corresponding
%   constant for each element of the input vector.
%
%   List of identifiers:
%       Generic astronomical constants:
%           1   Universal gravity constant (G) (from DITAN and Horizon) [km^3/(kg*s^2)]
%           2   Astronomical Unit (AU) (from DE405) [km]
%               Note:  The value for 1 au is from the IAU 2012 Resolution B1.
%       Sun related:
%           3   Sun mean radius (from DITAN) [km]
%           4   Sun planetary constant (mu = mass * G) (from DE405)
%               [km^3/s^2]
%           31  Energy flux density of the Sun (from Wertz,SMAD)
%               [W/m2 at 1 AU]
%       Other:
%           5   Speed of light in the vacuum (definition in the SI and Horizon) [km/s]
%           6   Standard free fall (the acceleration due to gravity on the
%               Earth's surface at sea level) (from Wertz,SMAD) [m/s^2]
%           7   Mean distance Earth-Moon (from Wertz,SMAD) [km]
%           8   Obliquity (angle) of the ecliptic at Epoch 2000 (from
%               Horizon) [rad]
%           9   Gravitatonal field constant of the Earth (from Wertz,SMAD,
%               taken from JGM-2). This should be used in conjunction to
%               Earth radius = 6378.1363 km
%           32  Days in a Julian year y = 365.25 d  (from Horizon)
%       Planetary constants of the planets (mu = mass * G) [km^3/s^2]:
%           11  Me      (from DE405)
%           12  V       (from DE405)
%           13  E       (from DE405)
%           14  Ma      (from DE405)
%           15  J       (from DE405)
%           16  S       (from DE405)
%           17  U       (from DE405)
%           18  N       (from DE405)
%           19  P       (from DE405)
%           20  Moon    (from DE405)
%       Mean radius of the planets [km]:
%           21  Me      (from Horizon)
%           22  V       (from Horizon)
%           23  E       (from Horizon)
%           24  Ma      (from Horizon)
%           25  J       (from Horizon)
%           26  S       (from Horizon)
%           27  U       (from Horizon)
%           28  N       (from Horizon)
%           29  P       (from Horizon)
%           30  Moon    (from Horizon)
%
%   Notes for upgrading this function:
%       It is possible to add new constants.
%       - DO NOT change the structure of the function, as well as its
%           prototype.
%       - DO NOT change the identifiers of the constants that have already
%           been defined in this function. If you want to add a new
%           constant, use an unused identifier.
%       - DO NOT add constants that can be easily computed starting form
%           other ones (avoid redundancy).
%       Contact the author for modifications.
%
% INPUT:
%   in      Vector of identifiers of required constants.
%
% OUTPUT:
%   out     Vector of constants.
%
% EXAMPLE:
%   astroConstants([2, 4, 26])
%      Returns a row vector in which there is the value of the AU, the Sun
%      planetary constant and the mean radius of Saturn.
%
%   astroConstants(10 + [1:9])
%      Returns a row vector with the planetary constant of each planet.
%
% REFERENCES:
%   - DITAN (Direct Interplanetary Trajectory Analysis), Massimiliano
%       Vasile, 2006.
%	- Wertz J. R., Larson W. J., "Space Mission Analysis and Design", Third
%       Edition, Space Technology Library 2003.
%   [A]   DE405 - http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf
%   [B]   Explanatory Supplement to the Astronomical Almanac. 1992. K. P.
%         Seidelmann, Ed., p.706 (Table 15.8) and p.316 (Table 5.8.1),
%         University Science Books, Mill Valley, California. 
%   [C]   Tholen, D.J. and Buie, M.W. 1990. "Further Analysis of
%         Pluto-Charon Mutual Event Observations" BAAS 22(3):1129.
%   [D]   Seidelmann, P.K. et al. 2007. "Report of the IAU/IAG Working
%         Group on cartographic coordinates and rotational elements: 2006"
%         Celestial Mech. Dyn. Astr. 98:155-180. 
%   [F]   Anderson, J.D., et al. 1987. "The mass, gravity field, and
%         ephemeris of Mercury" Icarus 71:337-349.
%   [G]   Konopliv, A.S., et al. 1999. "Venus gravity: 180th degree and
%         order model" Icarus 139:3-18.
%   [H]   Folkner, W.M. and Williams, J.G. 2008. "Mass parameters and
%         uncertainties in planetary ephemeris DE421." Interoffice Memo.
%         343R-08-004 (internal document), Jet Propulsion Laboratory,
%         Pasadena, CA. 
%   [I]   Jacobson, R.A. 2008. "Ephemerides of the Martian Satellites -
%         MAR080" Interoffice Memo. 343R-08-006 (internal document),
%         Jet Propulsion Laboratory, Pasadena, CA. 
%   [J]   Jacobson, R.A. 2005. "Jovian Satellite ephemeris - JUP230"
%         private communication. 
%   [K]   Jacobson, R.A., et al. 2006. "The gravity field of the Saturnian
%         system from satellite observations and spacecraft tracking data"
%         AJ 132(6):2520-2526. 
%   [L]   Jacobson, R.A. 2007. "The gravity field of the Uranian system and
%         the orbits of the Uranian satellites and rings" BAAS 39(3):453. 
%   [M]   Jacobson, R.A. 2008. "The orbits of the Neptunian satellites and
%         the orientation of the pole of Neptune" BAAS 40(2):296. 
%   [N]   Jacobson, R.A. 2007. "The orbits of the satellites of Pluto -
%         Ephemeris PLU017" private communication.
%   [W1]  http://ssd.jpl.nasa.gov/?planet_phys_par Last retrieved
%         20/03/2013
%   [W2]  http://ssd.jpl.nasa.gov/?sat_phys_par Last retrieved
%         20/03/2013
%   [W3]  http://ssd.jpl.nasa.gov/horizons.cgi Last retrieved
%         20/03/2013
%   [M1]  Bills, B.G. and Ferrari, A.J. 1977. ``A Harmonic Analysis of
%         Lunar Topography'', Icarus 31, 244-259.
%   [M2]  Standish, E. M. 1998. JPL Planetary and Lunar Ephemerides,
%         DE405/LE405.
%   [M3]  Lunar Constants and Models Document, Ralph B. Roncoli, 23 Sept 2005,
%         JPL Technical Document D-32296 
%
%
% CALLED FUNCTIONS:
%   (none)
%
% AUTHOR:
%   Matteo Ceriotti, 2006, MATLAB, astroConstants.m
%
% PREVIOUS VERSION:
%   Matteo Ceriotti, 2006, MATLAB, astro_constants.m, Ver. 1.2
%       - Header and function name in accordance with guidlines.
%
% CHANGELOG:
%   26/10/2006, Camilla Colombo: Updated.
%   22/10/2007, Camilla Colombo: astroConstants(8) added (Obliquity (angle)
%       of the ecliptic at Epoch 2000).
%   02/10/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%   12/11/2010, Camilla Colombo: astroConstants(9) added (J2) Note: the
%       present value of J2 is not consistent with the value of the Earth
%       radius. This value of J2 should be used in conjunction to Earth
%       radius = 6378.1363 km
%   19/03/2013, Camilla Colombo: constants updated to NASA JPL website.
%       References added.
%   20/03/2013, REVISION, Francesca Letizia.
%   22/03/2013, Francesca Letizia: all GM from DE405.
%
% -------------------------------------------------------------------------

% 9: J2
% 32: 365.25

out = zeros(1,length(in));
for i=1:length(in)
    switch in(i)
        case 1
            out(i)=6.67259e-20; % From DITAN and Horizon
        case 2
            out(i)=149597870.691; % From DE405
        case 3
            % out(i)=700000; % From DITAN
            out(i)=6.955*10^5; % From Horizon [W3]
        case 4
            % out(i)=0.19891000000000E+31*6.67259e-20; % From DITAN
            out(i)=1.32712440017987E+11; % From DE405 [A]
        case 5
            out(i)=299792.458; % Definition in the SI, Horizon, DE405
        case 6
            out(i)=9.80665; % Definition in Wertz, SMAD
        case 7
            % out(i)=384401; % Definition in Wertz, SMAD
            out(i)=384400; % From Horizon [W3]
        case 8
            % out(i)=23.43928111*pi/180; % Definition in Wertz, SMAD
            out(i)=84381.412/3600*pi/180; % Definition in Horizon
            % obliquity of ecliptic (J2000)    epsilon = 84381.412 (± 0.005) arcsec 
        case 9
            out(i)=0.1082626925638815e-2; % Definition in Wertz, SMAD
        case 11
            % out(i)=0.33020000000000E+24*6.67259e-20; % From DITAN
            %out(i)=0.330104E+24*6.67259e-20;    % From Horizon [F]
            out(i)=2.203208E+4;    % From DE405
        case 12
            % out(i)=0.48685000000000E+25*6.67259e-20; % From DITAN
            %out(i)=4.86732E+24*6.67259e-20;     % From Horizon [G]
            out(i)=3.24858599E+5; % From DE405
        case 13
            % out(i)=0.59736990612667E+25*6.67259e-20; % From DITAN
            % out(i)=5.97219E+24*6.67259e-20;     % From Horizon [H]
            out(i) = 3.98600433e+5; % From DE405
        case 14
            % out(i)=0.64184999247389E+24*6.67259e-20; % From DITAN
            %out(i)=0.641693E+24*6.67259e-20; 	% From Horizon [I]
            out(i) = 4.2828314E+4; %Frome DE405
        case 15
            % out(i)=0.18986000000000E+28*6.67259e-20; % From DITAN
            %out(i)=1898.13E+24*6.67259e-20; 	% From Horizon [J]
            out(i) = 1.26712767863E+08; % From DE405
        case 16
            % out(i)=0.56846000000000E+27*6.67259e-20; % From DITAN
            % out(i)=568.319E+24*6.67259e-20;     % From Horizon [k]
            out(i) = 3.79406260630E+07; % From DE405
        case 17
            % out(i)=0.86832000000000E+26*6.67259e-20; % From DITAN
            % out(i)=86.8103E+24*6.67259e-20;     % From Horizon [L]
            out(i)= 5.79454900700E+06; % From DE405
        case 18
            % out(i)=0.10243000000000E+27*6.67259e-20; % From DITAN
            % out(i)=102.410E+24*6.67259e-20;     % From Horizon [M]
            out(i) = 6.83653406400E+06; % From DE405
        case 19
            % out(i)=0.14120000000000E+23*6.67259e-20; % From DITAN
            %out(i)=.01309E+24*6.67259e-20;     % From Horizon [N]
            out(i) = 9.81601000000E+02; %From DE405
        case 20
            % out(i)=0.73476418263373E+23*6.67259e-20; % From DITAN
             out(i)=4902.801;                 % From Horizon  [M2]
            %out(i)=4902.801076;                % From Horizon  [M3]
        case 21
            % out(i)=0.24400000000000E+04; % From DITAN
            out(i)=2439.7; % From Horizon [D]
        case 22
            % out(i)=0.60518000000000E+04; % From DITAN
            out(i)=6051.8; % From Horizon [D]
        case 23
            % out(i)=0.63781600000000E+04; % From DITAN
            % out(i)=6371.00; % From Horizon [B]
            out(i)=6371.01; % From Horizon [W3]
        case 24
            % out(i)=0.33899200000000E+04; % From DITAN
            % out(i)=3389.50; % From Horizon [D]
            out(i)=3389.9; % From Horizon [W3]            
        case 25
            % out(i)=0.69911000000000E+05; % From DITAN
            out(i)=69911;   % From Horizon [D]
        case 26
            % out(i)=0.58232000000000E+05; % From DITAN
            out(i)=58232;   % From Horizon [D]
        case 27
            % out(i)=0.25362000000000E+05; % From DITAN
            out(i)=25362;   % From Horizon [D]
        case 28
            % out(i)=0.24624000000000E+05; % From DITAN
            % out(i)=24622;   % From Horizon [D]
            out(i)= 24624; % From Horizon [W3]            
        case 29
            % out(i)=0.11510000000000E+04; % From DITAN
            out(i)=1151; 	% From Horizon [C]
        case 30
            % out(i)=0.17380000000000E+04; % From DITAN
            % out(i)=1737.5;  % From Horizon [M1]
            out(i)=1738.0;    % From Horizon  [M3]
        case 31
            out(i)=1367; % From Wertz, SMAD
            % out(i)=1367.6;  % From Horizon  [W3]
        case 32
            out(i)=365.25; % From Horizon
        % Add an identifier and constant here. Prototype:
        % case $identifier$
        %     out(i)=$constant_value$;
        otherwise
            warning('Constant identifier %d is not defined!',in(i));
            out(i)=0;
    end
end

end