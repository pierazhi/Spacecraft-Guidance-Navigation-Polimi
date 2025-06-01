% Spacecraft Guidance & Navigation
% Assignment # 1, Exercise 1
% Author: Piercarlo Fontana

%% 1) Lagrangian Poinst
clc; clearvars; close all; format long g

% Define symbolic variables
syms x y z vx vy vz mu real

% Define the distance to the primary bodies
r1 = sqrt((x + mu)^2 + y^2 + z^2);
r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

% Define the potential energy function
U = 0.5*(x^2 + y^2) + (1 - mu)/r1 + mu/r2 + 0.5*mu*(1 - mu);

k = (1-mu)/r1^3+mu/r2^3;

% Compute the gradients (partial derivatives)
dUdx = x + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(3/2)) - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(3/2));
dUdy = y - (mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2);
dUdz = (z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) - (mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2);

% Substitute the value of mu (for the Earth-Moon system)
mu_value = 0.012150;

% Define the system of equations (setting the partial derivatives to zero)
eq1 = subs(dUdx, mu, mu_value) == 0;
eq2 = subs(dUdy, mu, mu_value) == 0;
eq3 = subs(dUdz, mu, mu_value) == 0;

% Solve the equations for equilibrium points
solutions = solve([eq1, eq2, eq3], [x, y, z]);

% Display the solutions
disp('Solutions for the Lagrange points:');
solutions_x = vpa(solutions.x, 11);
solutions_y = vpa(solutions.y, 11);
solutions_z = vpa(solutions.z, 11);

% Print the Lagrange points in the correct order
fprintf('Lagrange Points:\n');
fprintf('L1: x = %.10f, y = %.10f, z = %.10f\n', solutions_x(5), solutions_y(5), solutions_z(5)); % L1
fprintf('L2: x = %.10f, y = %.10f, z = %.10f\n', solutions_x(2), solutions_y(2), solutions_z(2)); % L2
fprintf('L3: x = %.10f, y = %.10f, z = %.10f\n', solutions_x(1), solutions_y(1), solutions_z(1)); % L3
fprintf('L4: x = %.10f, y = %.10f, z = %.10f\n', solutions_x(4), solutions_y(4), solutions_z(4)); % L4
fprintf('L5: x = %.10f, y = %.10f, z = %.10f\n\n', solutions_x(3), solutions_y(3), solutions_z(3)); % L5

% Compute the Jacobi Constant for each Lagrange point
C = zeros(1, 5);  % Initialize Jacobi Constant array
cont = 1;

for i = [5 2 1 4 3]
    temp = Jacobi([solutions_x(i), solutions_y(i), solutions_z(i), 0, 0, 0], mu_value);
    C(cont) = temp;
    cont = cont + 1;
end

% Display the Jacobi Constants for each Lagrange point
fprintf('Jacobi Constants:\n');
for i = 1:5
    fprintf('L%d: C = %.10f\n', i, C(i));
end
    
% Plot the Lagrange points in 3D
figure;
hold on;

% Earth and Moon positions
scatter(-mu_value, 0, 50, 'g', 'filled', 'DisplayName', 'Earth') % Earth at (-mu, 0)
scatter(1 - mu_value, 0, 50, 'b', 'filled', 'DisplayName', 'Moon') % Moon at (1 - mu, 0)

% Add floating labels
text(-mu_value + 0.1, 0.01, 'Earth', 'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold')
text(1 - mu_value, -0.1, 'Moon', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold')

%plotPlanet(3, [mu 0 0], gca, 1/200000);
%plotPlanet(11, [1-mu 0 0], gca, 1/50000);

% Plot Lagrange points
plot3(solutions_x, solutions_y, solutions_z, '.', 'MarkerSize', 10, 'DisplayName', 'Lagrange Points');

% Define the correct order of Lagrange points labels
lagrange_labels = {'L3', 'L2', 'L5', 'L4', 'L1'};

% Label Lagrange points with the new order
for i = 1:length(solutions_x)
    % Use the new lagrange_labels to label the points
    text(solutions_x(i), solutions_y(i), solutions_z(i), ...
        sprintf(' %s', lagrange_labels{i}), ...
        'FontSize', 15, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

% Set plot labels and title
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Lagrange Points @ EMRF');
grid on;
axis equal;

% Explicitly define the legend to ensure clarity
legend({'Earth', 'Moon', 'Lagrangian Points'}, 'Location', 'best');

% View the plot from a suitable angle
view(3);
hold off;

%% 2) - Periodic Halo Orbit with C = 3.09

clc; clearvars; format long g

% Initial Conditions
x0 = 1.068792441776;    
y0 = 0;
z0 = 0.071093328515;
vx0 = 0;
vy0 = 0.319422926485;
vz0 = 0;
mu = 0.012150; % Earth-Moon CRTBP 

% Calculate the radii and period of the orbit

r10 = sqrt((x0 + mu)^2 + y0^2 + z0^2);
r20 = sqrt((x0 + mu - 1)^2 + y0^2 + z0^2);
k = (1 - mu) / r10^3 + mu / r20^3;
T0 = 2 * pi / sqrt(k);

% Compute correction on the initial state using a pseudo-Newton method
err = 1;    % Initial error value
tol = 1e-14 ; % Desired tolerance
C_target = 3.09;
err_pos = 1;
err_Jac = 1;

% Initial guesses for the updated variables

vy0_new = vy0;
x0_new = x0;
z0_new = z0;
T = T0;

% Main correction loop
while err_pos > tol || err_Jac >tol 

    xx0_new=[x0_new;y0;z0_new;vx0;vy0_new;vz0];
    
    [xf,PHI,te] = propagate_3D(0,xx0_new,T,mu);
    Jf = Jacobi(xx0_new, mu);

    % Compute potential derivatives at xf
    r1 = sqrt((xf(1) + mu)^2 + xf(2)^2 + xf(3)^2);
    r2 = sqrt((xf(1) + mu - 1)^2 + xf(2)^2 + xf(3)^2);

    x = xf;

    dJdx = 2*x(1) + ((2*mu + 2*x(1))*(mu - 1))/((mu + x(1))^2 + x(2)^2 + x(3)^2)^(3/2) ...
       - (mu*(2*mu + 2*x(1) - 2))/((mu + x(1) - 1)^2 + x(2)^2 + x(3)^2)^(3/2);

    dJdz = (2*x(3)*(mu - 1))/((mu + x(1))^2 + x(2)^2 + x(3)^2)^(3/2) ...
       - (2*mu*x(3))/((mu + x(1) - 1)^2 + x(2)^2 + x(3)^2)^(3/2);
        
    dJdvy = -2*x(5);

    dUdx = xf(1) - (1-mu)/r1^3*(mu+xf(1)) + mu/r2^3*(1-mu-xf(1));
    dUdz= -(1-mu)*xf(3)/r1^3-mu*xf(3)/r2^3;
        
    % Compute the deviation in the final state ( the reference values are
    % equal to zero)

    vxf=xf(4);
    vzf=xf(6);

    vyf=xf(5);
    yf=xf(2);
    
    STM = [ PHI(2,1) , PHI(2,3), PHI(2,5) , vyf;
            PHI(4, 1), PHI(4,3), PHI(4, 5), dUdx + 2*vyf;
            PHI(6, 1), PHI(6,3), PHI(6, 5), dUdz;
            dJdx,      dJdz,     dJdvy,     0];  

    % Define the final states vector
    final_states = [-yf; -vxf; -vzf; Jf - C_target];

    % Update the state
    dxx = STM \ final_states;

    % Update variables
    x0_new = x0_new + dxx(1);
    z0_new = z0_new + dxx(2);
    vy0_new = vy0_new + dxx(3);
    T= te + dxx(4);
    
    err_pos=max(abs(vxf),abs(vzf));
    
    err_Jac = Jf - C_target;

end

xx0_new = [x0_new;y0;z0_new;vx0;vy0_new;vz0];

Jacobi(xx0_new,mu);

% Propagate the final corrected state of the halo orbit

[~, ~ , ~, xxc, ~] = propagate_3D(0, xx0_new, te, mu);

% Symmetry wrt y=0
xxc = symmetry(xxc);

fprintf('Corrected State for C = 3.09:\n');
fprintf('x0_new: %.10f\n', x0_new);
fprintf('y0_new: %.10f\n', y0);
fprintf('z0_new: %.10f\n', z0_new);
fprintf('vx0_new: %.10f\n', vx0);
fprintf('vy0_new: %.10f\n', vy0_new);
fprintf('vz0_new: %.10f\n', vz0);
fprintf('Final period (T): %.12f\n\n', T);

final_Jacobi = Jacobi(xx0_new, mu);
fprintf('Final Jacobi Constant: %.12f\n', final_Jacobi);

% Plot corrected trajectory
figure
hold on
plot3(xxc(:, 1), xxc(:, 2), xxc(:, 3), 'LineWidth', 2, 'Color', 'Magenta', 'DisplayName', 'Halo Orbit')
text(1.1556799131 + 0.01, 0.01, 'L2', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold')
scatter(1.1556799131, 0, 50, 'r', 'filled', 'DisplayName', 'L2');
xlabel('X [-]', 'FontSize', 14)
ylabel('Y [-]', 'FontSize', 14)
zlabel('Z [-]', 'FontSize', 14)
legend('show', 'Location', 'best', 'FontSize', 12)
title('Halo Orbit C = 3.09')
grid on
hold off
view(3)

% Plot corrected trajectory
figure
hold on
plot3(xxc(:, 1), xxc(:, 2), xxc(:, 3), 'LineWidth', 2, 'Color', 'Magenta', 'DisplayName', 'Halo Orbit')
text(1.1556799131 + 0.01, 0.01, 'L2', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold')
scatter(1.1556799131, 0, 50, 'r', 'filled', 'DisplayName', 'L2');
xlabel('x [-]', 'FontSize', 14)
ylabel('y [-]', 'FontSize', 14)
zlabel('z [-]', 'FontSize', 14)
legend('show', 'Location', 'best', 'FontSize', 12)
title('Halo Orbit C = 3.09')
grid on
hold off
view(2)

%% 3) - Families of Periodic Halo Orbits with C from 3.04 to 3.09
clc; clearvars; format long g

% Initial Conditions
x0 = 1.068792441776;    
y0 = 0;
z0 = 0.071093328515;
vx0 = 0;
vy0 = 0.319422926485;
vz0 = 0;
mu = 0.012150; % Earth-Moon CRTBP 

xx0 = [x0; y0; z0; vx0; vy0; vz0];

% Calculate the radii and period of the orbit
r10 = sqrt((x0 + mu)^2 + y0^2 + z0^2);
r20 = sqrt((x0 + mu - 1)^2 + y0^2 + z0^2);
k = (1 - mu) / r10^3 + mu / r20^3;
T0 = 2 * pi / sqrt(k);

tol = 1e-14; % Desired tolerance
Nmax = 100;  % Maximum number of iterations
err_pos = 1;
err_Jac = 1;
contatore = 1;
iter = 0;

% Initial guesses for the updated variables
vy0_new = vy0;
x0_new = x0;
z0_new = z0;

n_step = 10;
C_target = linspace(3.09, 3.04, n_step);

% Preallocate arrays
XFamily = zeros(length(C_target), 6);
C_final = zeros(length(C_target), 1);
T_family = zeros(length(C_target), 1);

for j = 1:length(C_target)

    fprintf('Starting calculation for C_target = %.6f\n', C_target(j));
    
    T = T0;
    err_pos = 1;
    iter = 0;
    err_Jac = 1;
    
    while err_pos > tol || err_Jac > tol 

            xx0_new=[x0_new;y0;z0_new;vx0;vy0_new;vz0];
    
            [xf,PHI,te] = propagate_3D(0,xx0_new,T,mu, true);
            Jf = Jacobi(xx0_new, mu);

            % Compute potential derivatives at xf
            r1 = sqrt((xf(1) + mu)^2 + xf(2)^2 + xf(3)^2);
            r2 = sqrt((xf(1) + mu - 1)^2 + xf(2)^2 + xf(3)^2);

            x = xf;

            dJdx = 2*x(1) + ((2*mu + 2*x(1))*(mu - 1))/((mu + x(1))^2 + x(2)^2 + x(3)^2)^(3/2) ...
                - (mu*(2*mu + 2*x(1) - 2))/((mu + x(1) - 1)^2 + x(2)^2 + x(3)^2)^(3/2);

            dJdz = (2*x(3)*(mu - 1))/((mu + x(1))^2 + x(2)^2 + x(3)^2)^(3/2) ...
                 - (2*mu*x(3))/((mu + x(1) - 1)^2 + x(2)^2 + x(3)^2)^(3/2);
        
            dJdvy = -2*x(5);

            dUdx = xf(1) - (1-mu)/r1^3*(mu+xf(1)) + mu/r2^3*(1-mu-xf(1));
            dUdz= -(1-mu)*xf(3)/r1^3-mu*xf(3)/r2^3;
        
             % Compute the deviation in the final state ( the reference values are
            % equal to zero)

            vxf=xf(4);
            vzf=xf(6);

            vyf=xf(5);
            yf=xf(2);
    
             STM = [ PHI(2,1) , PHI(2,3), PHI(2,5) , vyf;
            PHI(4, 1), PHI(4,3), PHI(4, 5), dUdx + 2*vyf;
            PHI(6, 1), PHI(6,3), PHI(6, 5), dUdz;
            dJdx,      dJdz,     dJdvy,     0];  

            % Define the final states vector
            final_states = [-yf; -vxf; -vzf; Jf - C_target(j)];

            % Update the state
            dxx = STM \ final_states;

            % Update variables
            x0_new = x0_new + dxx(1);
            z0_new = z0_new + dxx(2);
            vy0_new = vy0_new + dxx(3);
            T = te + dxx(4);
    
            err_pos=max(abs(vxf),abs(vzf));
    
            err_Jac = Jf - C_target(j);

            % Update iteration counter
                iter = iter+1;
            % compute the error 
    end

    sprintf('--Corrected State for C = %.3f\n--', C_target(j));
    fprintf('x0_new: %.10f\n', xx0_new(1));
    fprintf('y0_new: %.10f\n', xx0_new(2));
    fprintf('z0_new: %.10f\n', xx0_new(3));
    fprintf('vx0_new: %.10f\n', xx0_new(4));
    fprintf('vy0_new: %.10f\n', xx0_new(5));
    fprintf('vz0_new: %.10f\n\n', xx0_new(6));

    % Store the converged solution for this C_target
    XFamily(j,:) = [x0_new;y0;z0_new;vx0;vy0_new;vz0];
    T_family(j) = T;
end

for i = 1:length(C_target)
    C_final(i) = Jacobi(XFamily(i, :), mu); % Final Jacobi constant
end

figure
hold on
grid on

% Define a heatmap colormap (e.g., 'hot', 'cool', 'parula', etc.)
colormap('jet');  % Or use 'cool', 'parula', etc.
colors = colormap;  % Get the color map values
num_colors = size(colors, 1);  % Number of colors in the colormap

% Map Jacobi constants to the range of the colormap
C_min = min(C_target); % 3.04
C_max = max(C_target); % 3.09

% Create color indices based on the Jacobi constants
color_indices = round((C_target - C_min) / (C_max - C_min) * (num_colors - 1)) + 1;

for i = 1:size(XFamily, 1)
    % Propagate each orbit
    [~, ~, te0, xxc, ~] = propagate_3D(0, XFamily(i,:)', T_family(i), mu, true);
    
    % Symmetry wrt y=0
    xxc = symmetry(xxc);
    
    % Get the corresponding color from the colormap for this Jacobi constant
    plot_color = colors(color_indices(i), :);
    
    % Plot the orbit with the chosen color
    plot3(xxc(:, 1), xxc(:, 2), xxc(:, 3), 'Color', plot_color, 'LineWidth', 3, ...
          'DisplayName', sprintf('C = %.3f', C_target(i)));

end

% Mark L2 point
text(1.1556799131 + 0.01, 0.01, 'L2', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold')
scatter(1.1556799131, 0, 50, 'r', 'filled', 'DisplayName', 'L2');

% Customize the plot
xlabel('X [-]', 'FontSize', 14)
ylabel('Y [-]', 'FontSize', 14)
zlabel('Z [-]', 'FontSize', 14)
title('Family of Halo Orbits', 'FontSize', 14)
grid on
colorbar; % Add colorbar to show the gradient of Jacobi constants
clim([C_min C_max]) % Set colorbar limits to the range of Jacobi constants
view(3)

figure
hold on
grid on

% Define a heatmap colormap (e.g., 'hot', 'cool', 'parula', etc.)
colormap('jet');  % Or use 'cool', 'parula', etc.
colors = colormap;  % Get the color map values
num_colors = size(colors, 1);  % Number of colors in the colormap

% Map Jacobi constants to the range of the colormap
C_min = min(C_target); % 3.04
C_max = max(C_target); % 3.09

% Create color indices based on the Jacobi constants
color_indices = round((C_target - C_min) / (C_max - C_min) * (num_colors - 1)) + 1;

for i = 1:size(XFamily, 1)
    % Propagate each orbit
    [~, ~, te0, xxc, ~] = propagate_3D(0, XFamily(i,:)', T_family(i), mu, true);
    
    % Symmetry wrt y=0
    xxc = symmetry(xxc);
    
    % Get the corresponding color from the colormap for this Jacobi constant
    plot_color = colors(color_indices(i), :);
    
    % Plot the orbit with the chosen color
    plot3(xxc(:, 1), xxc(:, 2), xxc(:, 3), 'Color', plot_color, 'LineWidth', 3, ...
          'DisplayName', sprintf('C = %.3f', C_target(i)));
end

% Mark L2 point
%text(1.1556799131 + 0.01, 0.01, 'L2', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold')
%scatter(1.1556799131, 0, 50, 'r', 'filled', 'DisplayName', 'L2');

% Customize the plot
xlabel('x [-]', 'FontSize', 14)
ylabel('y [-]', 'FontSize', 14)
zlabel('z [-]', 'FontSize', 14)
title('Family of Halo Orbits', 'FontSize', 14)
grid on
colorbar; % Add colorbar to show the gradient of Jacobi constants
clim([C_min C_max]) % Set colorbar limits to the range of Jacobi constants
view(2)

% Save Figures

%saveAllFiguresToPath('/Users/pierazhi/Documents/MATLAB/GNC/Assignment 1 - Ex1/Images')

%% Functions

function C = Jacobi(xx, mu)
%--------------------------------------------------------------------------
% Function: Jacobi
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Computes the Jacobi constant for a given state vector in the 3D 
% Circular Restricted Three-Body Problem (CRTBP).
%
% The Jacobi constant is an energy-like conserved quantity used to analyze
% motion in the CR3BP and is derived from the pseudo-potential and kinetic energy.
%
% Inputs:
%   xx - State vector [x; y; z; vx; vy; vz] (6x1)
%   mu - Mass ratio of the two primary bodies
%
% Output:
%   C  - Jacobi constant (scalar)
%--------------------------------------------------------------------------

% Extract position and velocity components from the state vector
x  = xx(1);
y  = xx(2);
z  = xx(3);
vx = xx(4);
vy = xx(5);
vz = xx(6);

% Compute the pseudo-potential function Omega
Omega = 0.5*(x^2 + y^2) ...
      + (1 - mu)/sqrt((x + mu)^2 + y^2 + z^2) ...
      + mu/sqrt((x - 1 + mu)^2 + y^2 + z^2) ...
      + 0.5 * mu * (1 - mu);

% Calculate the Jacobi constant
C = 2*Omega - vx^2 - vy^2 - vz^2;



end

function [xf, PHIf, tf, xx, tt] = propagate_3D(t0, x0, tf, mu, varargin)
%--------------------------------------------------------------------------
% Function: propagate_3D
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Propagates the 3D state and State Transition Matrix (STM) of a body in the 
% Circular Restricted Three-Body Problem (CR3BP), with optional event detection.
%
% Integrates the equations of motion and STM using high-precision ODE integration.
% Useful for sensitivity analysis, trajectory optimization, and detection of
% symmetry crossings via event termination.
%
% Inputs:
%   t0       - Initial time
%   x0       - Initial state vector [x; y; z; vx; vy; vz] (6x1)
%   tf       - Final time
%   mu       - Mass ratio of the two primary bodies
%   varargin - Optional event flag (true to enable x-axis event detection)
%
% Outputs:
%   xf       - Final state vector (6x1)
%   PHIf     - Final STM (6x6 matrix)
%   tf       - Final time (possibly earlier if event triggered)
%   xx       - State history matrix (nx6)
%   tt       - Time vector (nx1)
%--------------------------------------------------------------------------

% Check if an event flag (evtFlag) is provided; default to true if not
if nargin > 4
    evtFlag = varargin{1};
else
    evtFlag = true;
end

% Calculate the time of flight (tof) as the difference between final and initial times
tof = tf - t0;

% Ensure the time of flight is positive
if tof <= 0
    error('Time of flight (tf - t0) must be greater than zero.');
end

% Initialize the State Transition Matrix (STM) as an identity matrix at t0
Phi0 = eye(6);

% Concatenate the initial state vector and STM for integration with ODE solver
x0Phi0 = [x0; Phi0(:)];

% Check if the initial state vector is the correct size
if numel(x0) ~= 6
    error('Initial conditions vector x0 must be a 6x1 vector.');
end

% Set options for the ODE solver with high precision and max step size,
% and configure event detection using x_axis_crossing function
options_STM = odeset('reltol', 1e-12, 'abstol', 1e-12, 'MaxStep', 0.01, ...
                     'Events', @(x,y) x_axis_crossing(x,y,evtFlag));

% Solve the equations of motion including the STM using ode78
[tt, xx] = ode113(@(t,x) xyCR3BP3D_STM(t,x,mu), [0 tof], x0Phi0, options_STM);

% Extract the final state vector, final STM, and final time
xf = xx(end, 1:6)';
PHIf = reshape(xx(end, 7:end), 6, 6);
tf = tt(end);

end

function [xx_sym] = symmetry(xx)
%--------------------------------------------------------------------------
% Function: symmetry
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Mirrors a semi-periodic orbit across the X-Z plane to generate the full
% periodic trajectory in the Circular Restricted Three-Body Problem (CRTBP).
%
% This function is useful when constructing periodic orbits using symmetry,
% such as Halo or Lyapunov trajectories, starting from a single half-period.
%
% Inputs:
%   xx      - Semi-period state matrix (l x 6), where each row is [x, y, z, vx, vy, vz]
%
% Output:
%   xx_sym  - Full-period state matrix (2l x 6), symmetric about X-Z plane
%--------------------------------------------------------------------------

% Determine the length of the input state vector (semi-period)
l = length(xx);
l2 = l * 2;

% Initialize the output matrix for the full-period state vector
xx_sym = zeros(2 * l, 6);

% Copy the first half (original semi-period state vector)
xx_sym(1:l, 1:6) = xx(:, 1:6);

% Mirror the state vector over the X-Z plane to complete the orbit
for i = 1:l
    xx_sym(l2 + 1 - i, 1) = xx(i, 1);   % x component remains the same
    xx_sym(l2 + 1 - i, 2) = -xx(i, 2);  % y component is mirrored
    xx_sym(l2 + 1 - i, 3) = xx(i, 3);   % z component remains the same
    xx_sym(l2 + 1 - i, 4) = -xx(i, 4);  % vx component is mirrored
    xx_sym(l2 + 1 - i, 5) = xx(i, 5);   % vy component remains the same
    xx_sym(l2 + 1 - i, 6) = -xx(i, 6);  % vz component is mirrored
end

end

function [value, isterminal, direction] = x_axis_crossing(~, xx, isTerminal)
%--------------------------------------------------------------------------
% Function: x_axis_crossing
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Event function to detect x-axis crossings (y = 0) in 3D CR3BP dynamics.
%
% Intended for use with MATLAB ODE solvers to enable event-based stopping
% or logging of trajectory crossings through the x-axis (i.e., y = 0).
%
% Inputs:
%   ~          - Unused time input
%   xx         - State vector [x; y; z; vx; vy; vz]
%   isTerminal - Boolean flag: stop integration on crossing if true
%
% Outputs:
%   value      - Event value (y component)
%   isterminal - Flag to stop integration when event is triggered
%   direction  - Direction of zero-crossing (0 = any)
%--------------------------------------------------------------------------
% Set value to y-coordinate for detecting x-axis crossings
value = xx(2);          % y component, which is zero at x-axis crossing
isterminal = isTerminal; % Stops integration if isTerminal is true
direction = 0;           % Detect crossing in both directions

end

function [dxdt] = xyCR3BP3D_STM(~, xx, mu)
%--------------------------------------------------------------------------
% Function: xyCR3BP3D_STM
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Computes the time derivative of the 3D CR3BP state and optionally the 
% State Transition Matrix (STM), depending on the length of the input vector.
%
% This function is used during the numerical propagation of both the state 
% and STM in the Circular Restricted Three-Body Problem (CR3BP). It forms 
% the right-hand side of the ODE integration used to compute the sensitivity 
% of the final state with respect to the initial conditions.
%
% Inputs:
%   ~    - Time variable (unused explicitly, but required by ODE solvers)
%   xx   - State vector (6x1) or state + STM vector (42x1: 6 + 6*6)
%   mu   - Mass ratio between the two primary bodies
%
% Output:
%   dxdt - Time derivative of the state (6x1) or of the full state + STM (42x1)
%--------------------------------------------------------------------------

 % Extract variables
    x  = xx(1);
    y  = xx(2);
    z = xx(3);
    vx = xx(4);
    vy = xx(5);
    vz = xx(6);

    % Put PHI in matrix form
    Phi = reshape(xx(7:end),6,6);
   
    % Compute distances from bodies 1 and 2
    r1 = sqrt((x + mu)^2 + y^2+z^2);
    r2 = sqrt((x + mu - 1)^2 + y^2+z^2);

   
    % Compute derivative of the potential
     dUdx = x - (1-mu)/r1^3*(mu+x) + mu/r2^3*(1-mu-x);
     dUdy = y - (1-mu)/r1^3*y - mu/r2^3*y;
     dUdz= -(1-mu)*z/r1^3-mu*z/r2^3;
     
    % Assemble the matrix dfdx 
    dfdx=[                                                                                                                                                                                                                              0,                                                                                                                                                                                     0,                                                                                                                                                                                 0,  1, 0, 0
                                                                                                                                                                                                                                        0,                                                                                                                                                                                     0,                                                                                                                                                                                 0,  0, 1, 0
                                                                                                                                                                                                                                        0,                                                                                                                                                                                     0,                                                                                                                                                                                 0,  0, 0, 1
    (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu + 2*x)*(mu + x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)) + 1,                                                                     (3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2),                                                                 (3*mu*z*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2),  0, 2, 0
                                                                                                        (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)), (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1,                                                                                   (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2), -2, 0, 0
                                                                                                        (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)),                                                                                       (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2), (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2),  0, 0, 0];
    % Compute the derivative of the STM
    Phidot = dfdx*Phi;

   

    % Assemble right-hand side
    dxdt = zeros(42,1);

    dxdt(1:3) = xx(4:6);
    dxdt(4)   =2*vy+dUdx;
    dxdt(5)   = -2*vx+dUdy;
    dxdt(6)    =dUdz;
    dxdt(7:end) = Phidot(:);
    
end