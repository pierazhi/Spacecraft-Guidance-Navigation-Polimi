% Spacecraft Guidance & Navigation
% Assignment # 1, Exercise 2
% Author: Piercarlo Fontana

%% 1) First guess solution
clearvars; clc; close all; set(0,'DefaultFigureVisible','on')
cspice_kclear
format long g

% Constants
const.mu = 1.21506683e-2;
const.Re = 6378e3;
const.hi = 167e3;
const.Rm = 1738e3;
const.hf = 100e3;
const.DU = 3.84405000e8;
const.ws = 0.925195985;
const.TU = 4.34256461;

const.ri = (const.Re + const.hi) / const.DU;
const.rf = (const.Rm + const.hf) / const.DU;

% Parameters
cspice_kclear
const.alpha = 0.2 * pi;
const.beta = 1.41;
const.delta = 4;
const.ti = 2;
const.tf = const.delta + const.ti;

const.r0 = (const.Re + const.hi) / const.DU;
const.v0 = const.beta * sqrt((1 - const.mu) / const.r0);

% Initial state
const.x0 = const.r0 * cos(const.alpha) - const.mu;
const.y0 = const.r0 * sin(const.alpha);
const.x0dot = -(const.v0 - const.r0) * sin(const.alpha);
const.y0dot = (const.v0 - const.r0) * cos(const.alpha);

fprintf('\n--- Initial Conditions ---\n');
fprintf('x0    = %12.6f DU\n', const.x0);
fprintf('y0    = %12.6f DU\n', const.y0);
fprintf('x0dot = %12.6f VU\n', const.x0dot);
fprintf('y0dot = %12.6f VU\n\n', const.y0dot);

const.X0 = [const.x0 const.y0 const.x0dot const.y0dot const.ti const.ti + const.delta]';

[~,~,~, xx, tt]  = propagate_PCRTBP(const.ti, const.X0(1:4), const.ti+const.delta, const.mu);   

XX_ECI = rotating2earthinertial(xx, tt, const.mu);

dpurple = [79, 121, 66]/norm([79, 121, 66]);

% Plot in the rotating frame
figure
hold on
grid on
axis equal
xlabel('X')
ylabel('Y')
title('Trajectory @ EMRF')
xlim([-1.5 3]);
plot(xx(:, 1), xx(:, 2), 'Color', dpurple, 'LineWidth', 2) % Main trajectory line

% Earth and Moon positions
scatter(-const.mu, 0, 50, dpurple, 'filled', 'DisplayName', 'Earth') % Earth at (-mu, 0)
scatter(1 - const.mu, 0, 50, 'b', 'filled', 'DisplayName', 'Moon') % Moon at (1 - mu, 0)

% Add floating labels
text(-const.mu + 0.1, 0.01, 'Earth', 'Color', dpurple, 'FontSize', 12, 'FontWeight', 'bold')
text(1 - const.mu + 0.1, 0.01, 'Moon', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold')

% Add legend
legend('Trajectory', 'Earth', 'Moon', 'Location', 'best')
hold off

% Plot in the Earth-centered frame
figure
hold on
grid on
axis equal
xlabel('X')
ylabel('Y')
title('Trajectory @ ECI')

% Draw Moon's orbit circle for reference
theta = linspace(0, 2 * pi, 100);
circle = (1 - const.mu) * [cos(theta); sin(theta)];
plot(circle(1, :), circle(2, :), 'LineWidth', 1, 'Color', 'blue', 'DisplayName', 'Moon Orbit') % Dashed line for orbit

% Plot trajectory
plot(XX_ECI(:, 1), XX_ECI(:, 2), 'Color', dpurple, 'LineWidth', 2, 'DisplayName', 'Trajectory')

% Earth as a reference point
scatter(0, 0, 80, dpurple, 'filled', 'DisplayName', 'Earth') % Earth at origin

% Add floating label for Earth
text(0.1, 0.1, 'Earth', 'Color', dpurple, 'FontSize', 12, 'FontWeight', 'bold')
% Add legend
legend('Moon Orbit', 'Trajectory', 'Earth', 'Location', 'best')
hold off

%% 2a) without providing derivative to the solver
clearvars -except const dpurple
format long g

LB(1:4) = -inf;
LB(5:6) = 0;
UB(1:4) = +inf;
UB(5) = 2*pi/const.ws;
UB(6) = UB(5) + 23*const.TU;

tic
options = optimset('Display', 'iter', 'LargeScale', 'off', 'Algorithm', 'active-set');
[Xsol, DV_NON_OPT] = fmincon(@(X) gradDV(X, const), const.X0, [], [], [], [], LB, UB, @(X) con_Grad(X, const), options);
t_esecuzione_NO_OPT = toc

DV_NON_OPT;
xx0_op = Xsol(1:4);
ti_op = Xsol(5);
tf_op = Xsol(6);

[~,~,~, xx, tt]  = propagate_PCRTBP(ti_op,xx0_op,tf_op,const.mu);   

fprintf('\n--- Without Derivatives ---\n');
fprintf('x0    = %12.6f DU\n', Xsol(1));
fprintf('y0    = %12.6f DU\n', Xsol(2));
fprintf('x0dot = %12.6f VU\n', Xsol(3));
fprintf('y0dot = %12.6f VU\n', Xsol(4));
fprintf('ti    = %12.3f TU\n', Xsol(5));
fprintf('tf    = %12.3f TU\n', Xsol(6));
fprintf('DV    = %12.4f km/s\n', DV_NON_OPT);
fprintf('TE    = %12.4f s\n\n', t_esecuzione_NO_OPT);

XX_ECI = rotating2earthinertial(xx, tt, const.mu);

const.sol_noder = XX_ECI;

dpurple = [79, 121, 66]/norm([79, 121, 66]);

% Plot in the rotating frame
figure
hold on
grid on
axis equal
xlabel('X')
ylabel('Y')
title('Trajectory @ EMRF')
plot(xx(:, 1), xx(:, 2), 'Color', dpurple, 'LineWidth', 2) % Main trajectory line

% Earth and Moon positions
scatter(-const.mu, 0, 50, dpurple, 'filled', 'DisplayName', 'Earth') % Earth at (-mu, 0)
scatter(1 - const.mu, 0, 50, 'b', 'filled', 'DisplayName', 'Moon') % Moon at (1 - mu, 0)

% Add floating labels
text(-const.mu + 0.1, 0.01, 'Earth', 'Color', dpurple, 'FontSize', 12, 'FontWeight', 'bold')
text(1 - const.mu + 0.1, 0.01, 'Moon', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold')

% Add legend
legend('Trajectory', 'Earth', 'Moon', 'Location', 'best')
hold off

% Plot in the Earth-centered frame
figure
hold on
grid on
axis equal
xlabel('X')
ylabel('Y')
title('Trajectory @ ECI')

% Draw Moon's orbit circle for reference
theta = linspace(0, 2 * pi, 100);
circle = (1 - const.mu) * [cos(theta); sin(theta)];
plot(circle(1, :), circle(2, :), 'LineWidth', 1, 'Color', 'blue', 'DisplayName', 'Moon Orbit') % Dashed line for orbit

% Plot trajectory
plot(XX_ECI(:, 1), XX_ECI(:, 2), 'Color', dpurple, 'LineWidth', 2, 'DisplayName', 'Trajectory')

% Earth as a reference point
scatter(0, 0, 80, dpurple, 'filled', 'DisplayName', 'Earth') % Earth at origin

% Add floating label for Earth
text(0.1, 0.1, 'Earth', 'Color', dpurple, 'FontSize', 12, 'FontWeight', 'bold')

% Add legend
legend('Moon Orbit', 'Trajectory', 'Earth', 'Location', 'best')
hold off

%% 2b) Providing derivative to the solver
clearvars -except const dpurple
format long g

LB(1:4) = -inf;
LB(5:6) = 0;
UB(1:4) = +inf;
UB(5) = 2*pi/const.ws;
UB(6) = UB(5) + 23*const.TU;

options_grad = optimoptions(@fmincon, 'Display', 'iter','Algorithm', 'active-set', ...
    'EnableFeasibilityMode', true,'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true);
tic
[Xsol, DV] = fmincon(@(X) gradDV(X, const), const.X0, [], [], [], [], LB, UB, @(X) con_Grad(X, const), options_grad);
t_OPT = toc;

DV;
xx0_op = Xsol(1:4);
ti_op = Xsol(5);
tf_op = Xsol(6);

fprintf('\n--- With Derivatives ---\n');
fprintf('x0    = %12.6f DU\n', Xsol(1));
fprintf('y0    = %12.6f DU\n', Xsol(2));
fprintf('x0dot = %12.6f VU\n', Xsol(3));
fprintf('y0dot = %12.6f VU\n', Xsol(4));
fprintf('ti    = %12.3f TU\n', Xsol(5));
fprintf('tf    = %12.3f TU\n', Xsol(6));
fprintf('DV    = %12.4f km/s\n', DV);
fprintf('TE    = %12.4f s\n\n', t_OPT);

[~,~,~, xx, tt]  = propagate_PCRTBP(ti_op,xx0_op,tf_op,const.mu);

XX_ECI = rotating2earthinertial(xx, tt, const.mu);

const.sol_der = XX_ECI;
const.SOL35 = [XX_ECI(1, 1:4)';ti_op;tf_op];

% Define custom colors
dpurple = [0.4660 0.6740 0.1880];

% Plot in the rotating frame
figure
hold on
grid on
axis equal
xlabel('X')
ylabel('Y')
title('Trajectory @ EMRF')
plot(xx(:, 1), xx(:, 2), 'Color', dpurple, 'LineWidth', 2) % Main trajectory line

% Earth and Moon positions
scatter(-const.mu, 0, 50, dpurple, 'filled', 'DisplayName', 'Earth') % Earth at (-mu, 0)
scatter(1 - const.mu, 0, 50, 'b', 'filled', 'DisplayName', 'Moon') % Moon at (1 - mu, 0)

% Add floating labels
text(-const.mu + 0.1, 0.01, 'Earth', 'Color', dpurple, 'FontSize', 12, 'FontWeight', 'bold')
text(1 - const.mu + 0.1, 0.01, 'Moon', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold')
% Add legend
legend('Trajectory', 'Earth', 'Moon', 'Location', 'best')
hold off

% Plot in the Earth-centered frame
figure
hold on
grid on
axis equal
xlabel('X')
ylabel('Y')
title('Trajectory @ ECI')

% Draw Moon's orbit circle for reference
theta = linspace(0, 2 * pi, 100);
circle = (1 - const.mu) * [cos(theta); sin(theta)];
plot(circle(1, :), circle(2, :), 'LineWidth', 1, 'Color', 'blue', 'DisplayName', 'Moon Orbit') % Dashed line for orbit

% Plot trajectory
plot(XX_ECI(:, 1), XX_ECI(:, 2), 'Color', dpurple, 'LineWidth', 2, 'DisplayName', 'Trajectory')

% Earth as a reference point
scatter(0, 0, 80, dpurple, 'filled', 'DisplayName', 'Earth') % Earth at origin

% Add floating label for Earth
text(0.1, 0.1, 'Earth', 'Color', dpurple, 'FontSize', 12, 'FontWeight', 'bold')

% Add legend
legend('Moon Orbit', 'Trajectory', 'Earth', 'Location', 'best')
hold off

%% 3) Multiple Shooting
clearvars -except const dpurple
format long g

N = 4;
    
time = zeros(1, N);
    for j = 1:4
        time(j) = const.ti + (const.tf - const.ti)*(j-1) / (N - 1);
    end

% Use propagate_PCRTBP to obtain the full state history and final state
[xx01, PHIf, tf, xx, ~] = propagate_PCRTBP(time(1), const.X0(1:4), time(2), const.mu);
[xx02, PHIf, tf, xx, ~] = propagate_PCRTBP(time(2), xx01, time(3), const.mu);
[xx03, PHIf, tf, xx, ~] = propagate_PCRTBP(time(3), xx02, time(4), const.mu);

X0 = [const.X0(1:4)' xx01' xx02' xx03' time(1) time(end)];
UB = zeros(1, 18);
LB = zeros(1, 18);
UB(1:16) = Inf;
LB(1:16) = -Inf;
UB(17) = 2 * pi / const.ws;
UB(18) = inf;

tic
options_grad = optimoptions(@fmincon, 'Display', 'iter', 'Algorithm', 'active-set', 'SpecifyObjectiveGradient',  true, 'SpecifyConstraintGradient', true, 'MaxIterations', 5000, 'MaxFunctionEvaluations', 100000);
ysol = fmincon(@(X) gradDV_multiple(X), X0, [], [], [], [], LB, UB, @(X) con_Grad_multiple(X), options_grad);
tempo_esecuzione = toc;

DV_MULTIPLE = gradDV_multiple(ysol);
xx1_sol = ysol(1:4);
xx2_sol = ysol(5:8);
xx3_sol = ysol(9:12);
xx4_sol = ysol(13:16);
t_initial = ysol(17);
t_final = ysol(18);

fprintf('\n--- Multiple Shooting ---\n');
fprintf('x0    = %12.6f DU\n', xx1_sol(1));
fprintf('y0    = %12.6f DU\n', xx1_sol(2));
fprintf('x0dot = %12.6f VU\n', xx1_sol(3));
fprintf('y0dot = %12.6f VU\n', xx1_sol(4));
fprintf('ti    = %12.3f TU\n', t_initial);
fprintf('tf    = %12.3f TU\n', t_final);
fprintf('DV    = %12.3f km/s\n', DV_MULTIPLE);
fprintf('TE    = %12.3f s\n\n', tempo_esecuzione);

tt_sol = zeros(1,4);

for j = 1:4
        tt_sol(j) = ysol(17) + (ysol(18) - ysol(17))*(j-1) / (N - 1);
end

[~,~,~, xx1, tt1] = propagate_PCRTBP(tt_sol(1),xx1_sol',tt_sol(2),const.mu);
[~,~,~, xx2, tt2] = propagate_PCRTBP(tt_sol(2),xx2_sol',tt_sol(3),const.mu);
[~,~,~, xx3, tt3] = propagate_PCRTBP(tt_sol(3),xx3_sol',tt_sol(4),const.mu);

% Plot in the rotating frame
figure
hold on
grid on
axis equal
xlabel('X')
ylabel('Y')
title('Trajectory @ EMRF')

% Earth and Moon positions
scatter(-const.mu, 0, 50, 'g', 'filled', 'DisplayName', 'Earth') % Earth at (-mu, 0)
scatter(1 - const.mu, 0, 50, 'b', 'filled', 'DisplayName', 'Moon') % Moon at (1 - mu, 0)
plot(xx1(:, 1), xx1(:, 2), 'Color', 'g', 'LineWidth', 2) % Main trajectory line
plot(xx2(:, 1), xx2(:, 2), 'Color', 'k', 'LineWidth', 2) % Main trajectory line
plot(xx3(:, 1), xx3(:, 2), 'Color', 'r', 'LineWidth', 2) % Main trajectory line

% Earth as a reference point
scatter(0, 0, 80, 'g', 'filled', 'DisplayName', 'Earth') % Earth at origin

% Add floating label for Earth
text(0.1, 0.1, 'Earth', 'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold')

legend('Earth', 'Moon', 'Segment 1', 'Segment 2', 'Segment 3', 'Location', 'best')
hold off

xx_eci1 = rotating2earthinertial(xx1, tt1, const.mu);
xx_eci2 = rotating2earthinertial(xx2, tt2, const.mu);
xx_eci3 = rotating2earthinertial(xx3, tt3, const.mu);

const.sol_mult = [xx_eci1; xx_eci2; xx_eci3];
% Plot in the Earth-centered frame
figure
hold on
grid on
axis equal
xlabel('X')
ylabel('Y')
title('Trajectory @ ECI')

% Draw Moon's orbit circle for reference
theta = linspace(0, 2 * pi, 100);
circle = (1 - const.mu) * [cos(theta); sin(theta)];
plot(circle(1, :), circle(2, :), 'LineWidth', 1, 'Color', 'blue', 'DisplayName', 'Moon Orbit') % Dashed line for orbit

plot(xx_eci1(:, 1), xx_eci1(:, 2), 'Color', 'g', 'LineWidth', 2) % Main trajectory line
plot(xx_eci2(:, 1), xx_eci2(:, 2), 'Color', 'k', 'LineWidth', 2) % Main trajectory line
plot(xx_eci3(:, 1), xx_eci3(:, 2), 'Color', 'r', 'LineWidth', 2) % Main trajectory line

% Earth as a reference point
scatter(0, 0, 80, 'g' , 'filled', 'DisplayName', 'Earth') % Earth at origin

% Add floating label for Earth
text(0.1, 0.1, 'Earth', 'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold')

% Add legend for clarity
legend('Moon Orbit', 'Segment 1', 'Segment 2', 'Segment 3', 'Earth', 'Location', 'best')
hold off

%% 4) Nbody
clearvars -except const dpurple
format long g

cspice_furnsh('kernels/assignment02.tm')

frame = 'ECLIPJ2000'; 

x_initial = [const.SOL35(1:2)*3.84405000*10^5;0;const.SOL35(3:4)*1.02454018;0]; 
t_initial = const.SOL35(5);
t_final = const.SOL35(6);

% Define the starting epoch
et_start = cspice_str2et('2024-09-28T00:00:00'); % Ephemerides Time 

theta0 = wrapTo2Pi(-const.ws * t_initial);
ti_et = fzero(@(t) angdiff(find_epoch(t, frame), theta0), [et_start, et_start+24*3600*28]); 
ti_utc = cspice_et2utc(ti_et, 'C', 4);
tf_et = ti_et + (t_final - t_initial)*const.TU*60*60*24;
tf_utc = cspice_et2utc(tf_et, 'C', 4);
fprintf('Departure date: %s \n', ti_utc)
fprintf('Arrival date:   %s \n',tf_utc)

% 2024 OCT 12 21:14:06.1782

% Define list of celestial bodies:
labels = {'Sun';
          'Earth';
          'Moon'};

% select integration frame string (SPICE naming convention)
center = 'Earth';

% Initialize propagation data by writing this function
bodies = nbody_init(labels);

options = odeset('reltol', 1e-12, 'abstol', [1e-6*ones(1,3),1e-9*ones(1,3)]);
[tt, xx] = ode113(@(t,x) nbody_shift_rhs(t, x, bodies, frame, center), [ti_et tf_et], x_initial, options);

others = {'Earth','Moon'};
for i = 1:length(others)
    rr_planet = zeros(3, length(tt));
    for j = 1:length(tt)
        rr_planet(:,j) = cspice_spkpos(others{i},tt(j),frame,'NONE',center);
    end
    
end

figure
plot3(xx(:,1), xx(:,2), xx(:,3), 'LineWidth', 2, 'Color', dpurple);
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title (['Trajectory @ ',frame])
axis equal
hold on
scatter(0, 0, 80, dpurple , 'filled', 'DisplayName', 'Earth') % Earth at origin
text(0.3*10^5, 0.3*10^5, 'Earth', 'Color', dpurple, 'FontSize', 12, 'FontWeight', 'bold')
plot3(rr_planet(1,:), rr_planet(2,:), rr_planet(3,:), 'LineWidth', 2, 'Color', 'b');
hold off
grid on
legend('Trajectory', 'Earth', 'Moon Trajectory', 'Location', 'best')
view(2)

figure
plot3(xx(:,1), xx(:,2), xx(:,3), 'LineWidth', 2, 'Color', dpurple);
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title (['Trajectory @ ',frame])
axis equal
hold on
scatter(0, 0, 80, dpurple , 'filled', 'DisplayName', 'Earth') % Earth at origin
text(0.3*10^5, 0.3*10^5, 'Earth', 'Color', dpurple, 'FontSize', 12, 'FontWeight', 'bold')
plot3(rr_planet(1,:), rr_planet(2,:), rr_planet(3,:), 'LineWidth', 2, 'Color', 'b');
hold off
grid on
legend('Trajectory', 'Earth', 'Moon Trajectory', 'Location', 'best')
view(3)

fprintf('\n--- N-Body ---\n');
fprintf('x0    = %12.8f km\n', xx(1, 1));
fprintf('y0    = %12.8f km\n', xx(1, 2));
fprintf('x0dot = %12.8f km/s\n', xx(1, 4));
fprintf('y0dot = %12.8f km/s\n', xx(1, 5));

figure
hold on
plot3(xx(:,1), xx(:,2), xx(:,3), 'LineWidth', 2, 'Color', dpurple);
plot3(const.sol_noder(:, 1)*3.84405000*10^5, const.sol_noder(:, 2)*3.84405000*10^5, const.sol_noder(:, 3)*3.84405000*10^5, 'LineWidth', 2, 'Color', 'r')
plot3(const.sol_der(:, 1)*3.84405000*10^5, const.sol_der(:, 2)*3.84405000*10^5, const.sol_der(:, 3)*3.84405000*10^5, 'LineWidth', 2, 'Color', 'k')
plot3(const.sol_mult(:, 1)*3.84405000*10^5, const.sol_mult(:, 2)*3.84405000*10^5, const.sol_mult(:, 3)*3.84405000*10^5, 'LineWidth', 2, 'Color', 'm')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title (['Trajectories Comparison @ ',frame])
axis equal
scatter(0, 0, 80, dpurple , 'filled', 'DisplayName', 'Earth') % Earth at origin
text(0.3*10^5, 0.3*10^5, 'Earth', 'Color', dpurple, 'FontSize', 12, 'FontWeight', 'bold')
plot3(rr_planet(1,:), rr_planet(2,:), rr_planet(3,:), 'LineWidth', 2, 'Color', 'b');
hold off
grid on
legend('N-Body', 'No Derivative', 'Derivative', 'Multiple Shooting','Earth', 'Moon Trajectory', 'Location', 'best')
view(2)

% Save Figures

% saveAllFiguresToPath('/Users/pierazhi/Documents/MATLAB/GNC/Assignment 1 - Ex2/Images')

%% Functions

function [c, ceq, gradc, gradceq] = con_Grad(X, const)
%--------------------------------------------------------------------------
% Function: con_Grad
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Computes the constraint values and their gradients for trajectory optimization 
% in the planar Circular Restricted Three-Body Problem (CR3BP).
%
% This function defines the equality constraints for imposing circular orbit 
% conditions at departure and arrival, along with corresponding velocity 
% constraints. An inequality constraint ensures that the final time is 
% greater than the initial time. If requested, analytical gradients for use in 
% gradient-based optimization routines are also provided.
%
% Inputs:
%   X       - State vector [xi, yi, vxi, vyi, ti, tf]
%   const   - Struct with constants (mu: mass ratio, ri: initial radius, rf: final radius)
%
% Outputs:
%   c       - Inequality constraint vector (1x1): [ti - tf]
%   ceq     - Equality constraint vector (4x1): circular orbit and velocity matching
%   gradc   - Gradient of inequality constraint (6x1)
%   gradceq - Jacobian of equality constraints (6x4)
%--------------------------------------------------------------------------
    
    mu = const.mu; % Earth-Moon mass ratio (this value should be confirmed)
    
    % Calculate normalized distances
    ri = const.ri; % Normalized initial distance
    rf = const.rf; % Normalized final distance
    
    % Unpack state vector
    xi = X(1);  yi = X(2);    % Initial position
    vxi = X(3); vyi = X(4);   % Initial velocity
    ti = X(5);  tf = X(6);    % Initial and final times
    
    % Propagate state to final time
    [xxf, Phif, ~, ~] = propagate_PCRTBP(ti, X(1:4), tf, mu);
    
    % Unpack final state
    xf = xxf(1); yf = xxf(2);    % Final position
    vxf = xxf(3); vyf = xxf(4);  % Final velocity
    
    % Combine equality constraints
    C1 = (xi + mu)^2 + yi^2 - ri^2;                      % Initial circular orbit constraint
    C2 = (xi + mu)*(vxi - yi) + yi*(vyi + xi + mu);    % Initial velocity constraint
    C3 = (xf + mu - 1)^2 + yf^2 - rf^2;                 % Final circular orbit constraint
    C4 = (xf + mu - 1)*(vxf - yf) + yf*(vyf + xf + mu - 1);  % Final velocity constraint

    ceq = [C1; C2; C3; C4];

    % Inequality constraint: ensure final time > initial time
    c = ti - tf;
    
    if nargout > 2    
        % Gradient of inequality constraint
        gradc = [0; 0; 0; 0; 1; -1];  % 6x1 column vector
        
        % Gradients for C3
        dC3dx =  [2*mu + 2*xf - 2, 2*yf, 0, 0] * Phif;
        dC3dt1 = [2*mu + 2*xf - 2, 2*yf, 0, 0] * -Phif * xyPCRTBP_STM(ti, X, mu, 'state');
        dC3dt2 = [2*mu + 2*xf - 2, 2*yf, 0, 0] * xyPCRTBP_STM(tf, xxf, mu, 'state');
        
        % Gradients for C4
        dC4dx =  [vxf, vyf, (mu + xf - 1), yf] * Phif;
        dC4dt1 = [vxf, vyf, (mu + xf - 1), yf] * -Phif * xyPCRTBP_STM(ti, X, mu, 'state');
        dC4dt2 = [vxf, vyf, (mu + xf - 1), yf] * xyPCRTBP_STM(tf, xxf, mu, 'state');
        
        % Jacobian of equality constraints
        gradceq = [2*mu + 2*xi, 2*yi, 0, 0, 0, 0;
                   vxi, vyi, mu + xi, yi, 0, 0;
                   dC3dx, dC3dt1, dC3dt2;
                   dC4dx, dC4dt1, dC4dt2]';
    end
end

function [g, c, gradg, gradc] = con_Grad_multiple(X)
%--------------------------------------------------------------------------
% Function: con_Grad_multiple
% Author: Piercarlo Fontana
% Date: 12 November 2024
% 
% This function calculates the constraint values and their gradients for a 
% multiple-shooting trajectory optimization problem in the Earth-Moon 
% Circular Restricted Three-Body Problem (CR3BP) framework. It defines 
% both inequality and equality constraints, considering positions and 
% velocities for initial and final orbit heights, as well as trajectory 
% continuity requirements across multiple segments. The function also 
% provides gradient matrices for constraints, aiding in optimization 
% algorithms requiring analytical derivatives.
%--------------------------------------------------------------------------

     % System constants
    Re = 6378e3;     % Earth radius [m]
    hi = 167e3;      % Initial orbit height [m]
    DU = 3.84405000e8; % Distance unit [m]
    Rm = 1738e3;     % Moon radius [m]
    hf = 100e3;      % Final orbit height [m]
    mu = 1.21506683e-2; % Earth-Moon mass ratio
    N = 4;
    
    % Time variables
    ti = X(17);
    tf = X(18);

    time = zeros(1, N);

    for j = 1:N
        time(j) = ti + (tf - ti)*(j-1) / (N - 1);
    end
    
    % Calculate inequality constraints (distance from Earth/Moon)

    eta1_e = (Re/DU)^2 - (X(1)+mu)^2 - X(2)^2;
    eta1_m = (Rm/DU)^2 - (X(1)+mu-1)^2 - X(2)^2;
    eta2_e = (Re/DU)^2 - (X(5)+mu)^2 - X(6)^2;
    eta2_m = (Rm/DU)^2 - (X(5)+mu-1)^2 - X(6)^2;
    eta3_e = (Re/DU)^2 - (X(9)+mu)^2 - X(10)^2;
    eta3_m = (Rm/DU)^2 - (X(9)+mu-1)^2 - X(10)^2;
    eta4_e = (Re/DU)^2 - (X(13)+mu)^2 - X(14)^2;
    eta4_m = (Rm/DU)^2 - (X(13)+mu-1)^2 - X(14)^2;
    
    % Time constraint
    tau = ti - tf;
    
    % Inequality Constraints

    g = [eta1_e; eta1_m; eta2_e; eta2_m; eta3_e; eta3_m; eta4_e; eta4_m; tau];

    % Initial and final conditions
    psi1_1 = (X(1)+mu)^2 + X(2)^2 - ((Re+hi)/DU)^2;
    psi1_2 = (X(1)+mu)*(X(3)-X(2)) + X(2)*(X(4)+X(1)+mu);
    psiN_1 = (X(13)+mu-1)^2 + X(14)^2 - ((Rm+hf)/DU)^2;
    psiN_2 = (X(13)+mu-1)*(X(15)-X(14)) + X(14)*(X(16)+X(13)+mu-1);

    % Extract states for each segment
    xx1 = X(1:4);
    xx2 = X(5:8);
    xx3 = X(9:12);
    xx4 = X(13:16);
    
    % Propagate trajectories
    [xxf_1, Phit1t2] = propagate_PCRTBP(time(1), xx1', time(2), mu);
    [xxf_2, Phit2t3] = propagate_PCRTBP(time(2), xx2', time(3), mu);
    [xxf_3, Phit3t4] = propagate_PCRTBP(time(3), xx3', time(4), mu);
        
    % Calculate matching conditions
    zeta1 =  xxf_1 - xx2';
    zeta2 =  xxf_2 - xx3';
    zeta3 =  xxf_3 - xx4';
    
    %Equality Constraints

    c = [zeta1; zeta2; zeta3; psi1_1; psi1_2; psiN_1; psiN_2];
    
if nargout>2
     % Assemble gradg - Inequality Constraints Gradient

     % Jacobian matrices for inequality constraints
     S1 = [-2*(X(1)+mu) -2*X(2) 0 0;
          -2*(X(1)+mu-1) -2*X(2) 0 0];
     S2 = [-2*(X(5)+mu) -2*X(6) 0 0;
          -2*(X(5)+mu-1) -2*X(6) 0 0];
     S3 = [-2*(X(9)+mu) -2*X(10) 0 0;
           -2*(X(9)+mu-1) -2*X(10) 0 0];
     S4 = [-2*(X(13)+mu) -2*X(14) 0 0;
           -2*(X(13)+mu-1) -2*X(14) 0 0];

     St = [1 -1];
    
     gradg = [S1 zeros(2, 14);
               zeros(2, 4)  S2 zeros(2, 10);
               zeros(2, 8)  S3 zeros(2, 6);
               zeros(2, 12) S4 zeros(2, 2);
               zeros(1, 16) St]';
     
    % Assemble gradc - Equality Constraints Gradient

    Q11 = - ((N-1) / (N-1)) * Phit1t2 * xyPCRTBP_STM(time(1), xx1, mu, 'state') + ((N-1-1) / (N - 1)) * xyPCRTBP_STM(time(2), xxf_1', mu, 'state');
    Q1N = -((1-1)/(N-1))*Phit1t2*xyPCRTBP_STM(time(1), xx1, mu, 'state') + (1/(N-1))*xyPCRTBP_STM(time(2), xxf_1', mu, 'state');

    Q21 = - ((N-2) / (N-1)) * Phit2t3 * xyPCRTBP_STM(time(2), xx2, mu, 'state') + ((N-2-1) / (N - 1)) * xyPCRTBP_STM(time(3), xxf_2', mu, 'state');
    Q2N = -((2-1)/(N-1))*Phit2t3*xyPCRTBP_STM(time(2), xx2, mu, 'state') + (2/(N-1))*xyPCRTBP_STM(time(3), xxf_2', mu, 'state');  

    Q31 = - ((N-3) / (N-1)) * Phit3t4 * xyPCRTBP_STM(time(3), xx3, mu, 'state') + ((N-3-1) / (N - 1)) * xyPCRTBP_STM(time(4), xxf_3', mu, 'state');
    Q3N = -((3-1)/(N-1))*Phit3t4*xyPCRTBP_STM(time(3), xx3, mu, 'state') + (3/(N-1))*xyPCRTBP_STM(time(4), xxf_3', mu, 'state');   
    
    % Add boundary condition gradients

     R1 = [2*(X(1)+mu) 2*X(2) 0 0;
           X(3) X(4) X(1)+mu X(2)];
     RN = [2*(X(13)+mu-1) 2*X(14) 0 0;
          X(15) X(16) X(13)+mu-1 X(14)];
     

     gradc = [Phit1t2 -eye(4) zeros(4) zeros(4) Q11 Q1N;
              zeros(4) Phit2t3 -eye(4) zeros(4) Q21 Q2N;
              zeros(4) zeros(4) Phit3t4 -eye(4) Q31 Q3N;
              R1       zeros(2, 14);
              zeros(2, 12) RN zeros(2)]';

end
end

function theta_x = find_epoch(et, frame)
%--------------------------------------------------------------------------
% Function: find_epoch
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Computes the angular position of the Sun with respect to the rotating frame 
% defined by the Earth-Moon Barycenter (EMB), using SPICE ephemerides.
%
% The function constructs a Moon-centered inertial frame aligned with the Moon's 
% position and velocity vectors and computes the angle of the Sun in this rotated frame.
%
% Inputs:
%   et    - Ephemeris time (seconds past J2000)
%   frame - String specifying the inertial frame (e.g., 'J2000')
%
% Output:
%   theta_x - Angular position of the Sun in the rotated frame [rad]
%--------------------------------------------------------------------------

    % Sun and Moon positions relative to EMB at epoch et
    r_S = cspice_spkezr('SUN',  et, frame, 'NONE', 'EMB');
    r_M = cspice_spkezr('MOON', et, frame, 'NONE', 'EMB');

    x_rot = r_M(1:3) / norm(r_M(1:3));
    z_rot = cross(r_M(1:3), r_M(4:6)) / norm(cross(r_M(1:3), r_M(4:6)));
    y_rot = cross(z_rot, x_rot);

    R = [x_rot, y_rot, z_rot]'; % Combine into a rotation matrix

    % Transform Sun's position to the rotated frame
    rS_rot = R * r_S(1:3);

    % Compute the angular position of the Sun in the rotated frame
    xS_rot = rS_rot(1);
    yS_rot = rS_rot(2);

    theta_x = wrapTo2Pi(atan2(yS_rot, xS_rot));

end

function [DV, gradDV] = gradDV(X, const)
%--------------------------------------------------------------------------
% Function: gradDV
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Calculates the total delta-V required for a two-impulse transfer trajectory 
% in the planar Circular Restricted Three-Body Problem (CR3BP) and optionally 
% computes its gradient with respect to the initial state and times.
%
% The delta-V is computed from velocity changes at the initial and final 
% points. Gradients are returned using STM-based sensitivity propagation.
%
% Inputs:
%   X       - State vector [xi, yi, vxi, vyi, ti, tf]
%   const   - Struct with constants (mu, ri, rf)
%
% Outputs:
%   DV      - Total delta-V scalar
%   gradDV  - Gradient of delta-V (6x1 vector)
%--------------------------------------------------------------------------
    
    mu = const.mu; % Earth-Moon mass ratio (this value should be confirmed)
    
    % Calculate normalized distances
    ri = const.ri; % Normalized initial distance
    rf = const.rf; % Normalized final distance
     
    % Extract state variables
    xi = X(1);
    yi = X(2);
    vxi = X(3); 
    vyi = X(4);

    ti = X(5);
    tf = X(6);

    % Calculate delta-V components
    DVi = sqrt((vxi - yi)^2 + (vyi + xi + mu)^2) - sqrt((1 - mu) / ri);
    
    % Propagate the state from initial to final time
    [xxf, Phif, ~, ~] = propagate_PCRTBP(X(5), X(1:4), X(6), mu);  
        
    xf = xxf(1);
    yf = xxf(2);
    vxf = xxf(3);
    vyf = xxf(4);
    
    DVf = sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2) - sqrt(mu / rf);
    
    % Total delta-V
    DV = DVi + DVf;
    
    % Calculate the gradient with respect to initial conditions
    if nargout > 1
         % Initialize gradient
        gradDV = zeros(6, 1); % Ensure gradient is 6-by-1
        
        dVidxi = 1 / sqrt( (vxi - yi)^2 + (vyi + xi + mu)^2)* [vyi + xi + mu;
                                                     yi - vxi;
                                                     vxi - yi;
                                                     vyi + xi + mu];

        dVfdxf = 1 / sqrt((vxf - yf)^2 + (mu + vyf + xf - 1)^2)* [mu + vyf + xf - 1;
                                                                 -vxf + yf;
                                                                 vxf - yf;
                                                                  mu + vyf + xf - 1];
        % Gradient with respect to initial state
        gradDV(1:4) = dVidxi + Phif'*dVfdxf;
        
        % Incorporate State Transition Matrix effects
        gradDV(5) =  dVfdxf' * -Phif * xyPCRTBP_STM(ti, X, mu, 'state');
        gradDV(6) =  dVfdxf' * xyPCRTBP_STM(tf,xxf, mu, 'state');
    end
end

function [J, gradJ] = gradDV_multiple(X)
%--------------------------------------------------------------------------
% Function: gradDV_multiple
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Computes the total delta-V and its gradient for a multi-leg trajectory 
% in the planar Circular Restricted Three-Body Problem (CR3BP).
%
% The function evaluates the velocity changes at the first and last segments 
% of the trajectory and, if requested, returns the gradient of the total 
% delta-V with respect to the relevant state variables.
%
% Inputs:
%   X      - Optimization vector including state components of multiple segments
%
% Outputs:
%   J      - Total delta-V cost function (scalar)
%   gradJ  - Gradient of J with respect to X (18x1 vector)
%--------------------------------------------------------------------------

 % Define constants
    Re = 6378e3; % Radius of Earth in meters
    hi = 167e3;  % Initial altitude in meters
    DU = 3.84405000e8; % Distance Unit: Distance from Earth to Moon in meters
    Rm = 1738e3; % Radius of Moon in meters
    hf = 100e3;  % Final altitude in meters
    mu = 1.21506683e-2; % Earth-Moon mass ratio (this value should be confirmed)
    
    % Calculate normalized distances
    ri = (Re + hi) / DU; % Normalized initial distance
    rf = (Rm + hf) / DU; % Normalized final distance
  
    % Extract state variables
    xi = X(1);
    yi = X(2);
    vxi = X(3); 
    vyi = X(4);

    xn = X(13);
    yn = X(14);
    vxn = X(15); 
    vyn = X(16);
    
    % Calculate delta-V components
    DVi = sqrt((vxi - yi)^2 + (vyi + xi + mu)^2) - sqrt((1 - mu) / ri);
    DVf = sqrt((vxn - yn)^2 + (vyn + xn + mu - 1)^2) - sqrt(mu / rf);
    
    J =  DVi + DVf;
    
    if nargout>1
        gradJ = zeros(18, 1);
        
        gradJ(1:4) = 1 / sqrt( (vxi - yi)^2 + (vyi + xi + mu)^2)* [vyi + xi + mu;
                                                     yi - vxi;
                                                     vxi - yi;
                                                     vyi + xi + mu];

        gradJ(13:16) = 1/sqrt((vxn - yn)^2 + (vyn + mu + xn - 1)^2) * [vyn + xn + mu - 1;
                                                            yn - vxn;
                                                            vxn - yn;
                                                          vyn + xn + mu - 1];
    end

end

function [xf, PHIf, tf, xx, tt] = propagate_PCRTBP(ti, x0, tf, mu)
%--------------------------------------------------------------------------
% Function: propagate_PCRTBP
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Propagates the state and State Transition Matrix (STM) of a spacecraft 
% in the planar Circular Restricted Three-Body Problem (CR3BP) between 
% specified initial and final times.
%
% The STM is integrated along with the state for sensitivity analysis or 
% optimization tasks. Returns full trajectory data and final STM.
%
% Inputs:
%   ti   - Initial time
%   x0   - Initial state vector [xi, yi, vxi, vyi]
%   tf   - Final time
%   mu   - Mass ratio parameter
%
% Outputs:
%   xf   - Final state vector at time tf
%   PHIf - Final STM at time tf (4x4 matrix)
%   tf   - Final time (echoed back)
%   xx   - State history over time
%   tt   - Time vector
%--------------------------------------------------------------------------

    % Initialize State Transition Matrix at t0
    Phi0 = eye(4);

    % Combine initial conditions with STM for ODE integration
    x0Phi0 = [x0; Phi0(:)];

    options = odeset('AbsTol', 1e-13, 'RelTol', 1e-12);
    [tt, xx] = ode113(@(t, x) xyPCRTBP_STM(t, x, mu, 'full'), [ti tf], x0Phi0, options);

    % Extract final state vector and State Transition Matrix
    xf = xx(end, 1:4)';
    PHIf = reshape(xx(end, 5:end), 4, 4);
    tf = tt(end);
end

function xx_earth_inertial = rotating2earthinertial(xx_rotating, t, mu)
%--------------------------------------------------------------------------
% Function: rotating2earthinertial
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Converts state vectors from the rotating CR3BP frame to the Earth-centered 
% inertial frame, accounting for rotation at time t.
%
% This transformation is essential for visualizing or analyzing trajectories 
% in inertial coordinates and is based on standard CR3BP frame rotation.
%
% Inputs:
%   xx_rotating - State vector(s) in rotating frame [x, y, x_dot, y_dot] (nx4)
%   t           - Time vector corresponding to each state row (nx1)
%   mu          - Earth-Moon mass ratio
%
% Outputs:
%   xx_earth_inertial - State vector(s) in Earth-centered inertial frame (nx4)
%--------------------------------------------------------------------------

    % Extract the variables
    x_rotating = xx_rotating(:,1);
    y_rotating = xx_rotating(:,2);
    x_dot_rotating = xx_rotating(:,3);
    y_dot_rotating = xx_rotating(:,4);

    % Compute the state in the Earth-centered inertial frame
    x_earth_inertial = (x_rotating + mu) .* cos(t) - y_rotating .* sin(t);
    y_earth_inertial = (x_rotating + mu) .* sin(t) + y_rotating .* cos(t);
    x_dot_earth_inertial = (x_dot_rotating - y_rotating) .* cos(t) - (y_dot_rotating + x_rotating + mu) .* sin(t);
    y_dot_earth_inertial = (x_dot_rotating - y_rotating) .* sin(t) + (y_dot_rotating + x_rotating + mu) .* cos(t);

    xx_earth_inertial = [x_earth_inertial, y_earth_inertial, x_dot_earth_inertial, y_dot_earth_inertial];
end

function [dxdt] = xyPCRTBP_STM(t, xx, mu, returnMode)
%--------------------------------------------------------------------------
% XYPCRTBP_STM Computes the state dynamics and State Transition Matrix (STM) 
% in the planar Circular Restricted Three-Body Problem (CR3BP) including 
% solar gravitational effects.
%
% Author: Piercarlo Fontana
% Date: 12 November 2024
%
% Inputs:
%   t          - Time variable
%   xx         - State vector with optional STM (4+16 = 20x1 for 'full' mode)
%   mu         - Mass ratio (gravitational parameter)
%   returnMode - Mode selection, 'state' for state dynamics or 'full' for STM
%
% Outputs:
%   dxdt       - State derivative or full derivative with STM
%--------------------------------------------------------------------------   
        ms = 3.28900541 * 10^5;
        rho = 3.88811143 * 10^2;
        ws = -9.25195985*10^-1;

        % Extract variables
        x  = xx(1);
        y  = xx(2);
        vx = xx(3); 
        vy = xx(4);

         % Compute derivatives of the potential

         dUdx = x - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(3/2)) - (ms*cos(t*ws))/rho^2 + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (ms*(2*x - 2*rho*cos(t*ws)))/(2*((x - rho*cos(t*ws))^2 + (y - rho*sin(t*ws))^2)^(3/2));
         dUdy = y - (ms*sin(t*ws))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (ms*(2*y - 2*rho*sin(t*ws)))/(2*((x - rho*cos(t*ws))^2 + (y - rho*sin(t*ws))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);

        switch returnMode
            case 'state'

                dxdt = zeros(4, 1);
                dxdt(1:2) = xx(3:4);
                dxdt(3) = dUdx + 2*vy;
                dxdt(4) = dUdy - 2*vx;

        case 'full' 
                % Reshape PHI into 4x4 matrix
                Phi = reshape(xx(5:end), 4, 4);
    
                dUdxx = (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(ws*t))^2 + (y - rho*sin(ws*t))^2)^(3/2) + (3*ms*(2*x - 2*rho*cos(ws*t))^2)/(4*((x - rho*cos(ws*t))^2 + (y - rho*sin(ws*t))^2)^(5/2)) + (3*mu*(2*mu + 2*x - 2)^2)/(4*((mu + x - 1)^2 + y^2)^(5/2)) - (3*(2*mu + 2*x)^2*(mu - 1))/(4*((mu + x)^2 + y^2)^(5/2)) + 1;    
                dUdyy = (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(ws*t))^2 + (y - rho*sin(ws*t))^2)^(3/2) + (3*ms*(2*y - 2*rho*sin(ws*t))^2)/(4*((x - rho*cos(ws*t))^2 + (y - rho*sin(ws*t))^2)^(5/2)) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1;
                dUdxy = (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) + (3*ms*(2*x - 2*rho*cos(t*ws))*(2*y - 2*rho*sin(t*ws)))/(4*((x - rho*cos(t*ws))^2 + (y - rho*sin(t*ws))^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2));
       
                % Assemble the matrix A(t) = dfdx (4x4 matrix)
                dfdx = [
                        0,    0,  1,  0;
                        0,    0,  0,  1;
                        dUdxx, dUdxy, 0, 2;
                        dUdxy, dUdyy, -2, 0
                        ];
    
                % Compute the derivative of the STM

                Phidot = dfdx * Phi;
                dxdt = zeros(20, 1);
                dxdt(1:4) = [xx(3:4); dUdx + 2*vy; dUdy - 2*vx];
                dxdt(5:end) = Phidot(:);
                
            otherwise
                error('Invalid return mode. Use ''state'' or ''full''.');
        end

      
    
end

function [dxdt] = nbody_shift_rhs(t, x, bodies, frame, center)
%NBODY_RHS Evaluates the right-hand-side of a N-body propagator
%   Evaluates the right-hand-side of a newtonian N-body propagator.
%   The integration centre is the Solar-System-Barycentre (SSB) and only
%   Newtonian gravitational accelerations are considered.
%
% Inputs:
%   t      : [1,1] ephemeris time (ET SPICE), seconds past J2000 (TDB)
%   x      : [6,1] cartesian state vector wrt desired object
%   bodies : [1,n] cell-array created with function nbody_init
%
% Outputs:
%   dxdt   : [6,1] RHS, newtonian gravitational acceleration only
%

if not( strcmpi(frame, 'ECLIPJ2000') || strcmpi(frame, 'J2000') )
    msg = 'Invalid integration reference frame, select either J2000 or ECLIPJ2000';
    error(msg);
end

if not( any(strcmp(center,{cell2mat(bodies).name})))
    msg = 'Invalid center selected, select one of the bodies';
    error(msg);
end

% Initialize right-hand-side

dxdt = zeros(6,1);

% Position derivative is object's velocity

dxdt(1:3) = x(4:6);

% Extract the object position from state x

rr_b0_obj = x(1:3);

GM0 = bodies{strcmp(center, {cell2mat(bodies).name})}.GM;
dxdt(4:6) = -GM0*rr_b0_obj/norm(rr_b0_obj)^3;

for i = 1:length(bodies)

    if strcmp(bodies{i}.name, center)
        continue
    end


    % Retrieve position and velocity of i-th celestial body wrt desired
    % center in inertial frame:

    rv_b0_body = cspice_spkezr(bodies{i}.name, t, frame, 'NONE', center);   

    % Extract object position wrt. i-th celestial body

    rr_body_obj = rr_b0_obj - rv_b0_body(1:3);

    % Compute non-inertial terms as in the slides (d, rho, q, f):

    r = rr_b0_obj;
    d = rr_body_obj;
    rho = rv_b0_body(1:3);

    q = dot(r, r-2*rho)/dot(rho, rho);
    f = q*(3+3*q+q^2)/(1+(1+q)^1.5);

    aa_grav = -bodies{i}.GM * 1/norm(d)^3*(r+rho*f);

    % Sum up acceleration to right-hand-side

    dxdt(4:6) = dxdt(4:6) + aa_grav;

end

end

function [bodies] = nbody_init(labels)
%NBODY_INIT Initialize planetary data for n-body propagation
%   Given a set of labels of planets and/or barycentres, returns a
%   cell array populated with structures containing the body label and the
%   associated gravitational constant.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 26/09/2021
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
%   labels : [1,n] cell-array with object labels
%
% Outputs:
%   bodies : [1,n] cell-array with struct elements containing the following
%                  fields
%                  |
%                  |--bodies{i}.name -> body label
%                  |--bodies{i}.GM   -> gravitational constant [km**3/s**2]
%
%
% Prerequisites:
%   - MICE (Matlab SPICE)
%   - Populated kernel pool (PCK kernels)
%

% Initialize output
bodies = cell(size(labels));

% Loop over labels
for i = 1:length(labels)
    % Store body label
    bodies{i}.name = labels{i};
    % Store body gravitational constant
    bodies{i}.GM   = cspice_bodvrd(labels{i}, 'GM', 1);
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













