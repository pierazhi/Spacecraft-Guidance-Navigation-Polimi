% Spacecraft Guidance & Navigation
% Assignment # 2, Exercise 1
% Author: Piercarlo Fontana

%% 1) Plot the mean and the ellipses associated with LinCov and UT
clc; close all; clearvars; format longG

rng(42)

const.mu = 1.21506683e-2;
const.xi =  [-0.011965533749906 -0.017025663128129];
const.vi = [10.718855256727338 0.116502348513671];
const.X0 = [const.xi const.vi]';
const.ti = 1.282800225339865;
const.tf = 9.595124551366348;
const.P0 = [1.041e-15, 6.026e-17, 5.647e-16, 4.577e-15; ...
            6.026e-17, 4.287e-18, 4.312e-17, 1.855e-16; ...
            5.647e-16 4.312e-17 4.432e-16 1.455e-15; ...
            4.577e-15 1.855e-16 1.455e-15 2.822e-14];

dt = 5;
const.conf = 0.9973;
tspan = linspace(const.ti, const.tf, dt);
Phif = zeros(4,4, length(tspan));
xf = zeros(4, length(tspan));

P_LinCov(:, :, 1) = const.P0;
mean_LinCov(:, 1) = const.X0;

P_r_LinCov = zeros(2, 2, length(tspan)-1);
P_v_LinCov = zeros(2, 2, length(tspan)-1);
P_r_LinCov(:, :, 1) = const.P0(1:2, 1:2);
P_v_LinCov(:, :, 1) = const.P0(3:4, 3:4);

eigen_LinCov_r = zeros(1, dt-1);
eigen_LinCov_v = zeros(1, dt-1);
eigen_LinCov_r(1) = 3*sqrt(max(eig(P_r_LinCov(:, :, 1))));
eigen_LinCov_v(1) = 3*sqrt(max(eig(P_v_LinCov(:, :, 1))));

P_LinCov(:, :, 1) = const.P0;
for j = 2:dt
    [xf(:, j), Phif(:, :, j),~, xx, tt]  = propagate_PCRTBP(const.ti, [const.xi const.vi]', tspan(j), const.mu); 
    P_LinCov(:, :, j) = Phif(:, :, j)*const.P0*Phif(:, :, j)';
    %P_LinCov(:, :, j) = Phif(:, :, j)*P_LinCov(:, :, j-1)*Phif(:, :, j)';

    mean_LinCov(:, j) = xf(:, j);

    % Extract position and velocity covariance sub-matrices
    P_r_LinCov(:, :, j) = P_LinCov(1:2, 1:2, j);  % Position covariance
    P_v_LinCov(:, :, j) = P_LinCov(3:4, 3:4, j);  % Velocity covariance

    eigen_LinCov_r(j) = 3*sqrt(max(eig(P_r_LinCov(:, :, j))));
    eigen_LinCov_v(j) = 3*sqrt(max(eig(P_v_LinCov(:, :, j))));
    
end

% Unscented Transform

clc; close all; clearvars -except const mean_LinCov P_r_LinCov eigen_LinCov_r eigen_LinCov_v  P_v_LinCov

dt = 5;
tspan = linspace(const.ti, const.tf, dt);
Phif = zeros(4,4, length(tspan)-1);
P_UT(:, :, 1) = const.P0;
mean_UT(:, 1) = const.X0;

n = 4;  % State dimension
alpha = 1;
beta = 2;

c = alpha^2 * n;
W0_m = 1 - n / (alpha^2 * n);
W0_c = (2 - alpha^2 + beta) - n / (alpha^2*n);
Wi = (1 / (2*alpha^2*n))*ones(1, 2*n);

W_m = [W0_m Wi];
W_c = [W0_c Wi];
  
sigma_points(:,1) = const.X0;  % First sigma point is the mean state
dX0 = sqrtm(c * const.P0);

 for i = 1:n
        sigma_points(:, i+1) = const.X0 + dX0(:, i);  
        sigma_points(:, i+n+1) = const.X0 - dX0(:, i);  
 end

% Initialize xf to store the propagated states of all sigma points at each time step
xf = zeros(n, 2*n+1, length(tspan)-1);

P_r_UT = zeros(2, 2, length(tspan)-1);
P_v_UT = zeros(2, 2, length(tspan)-1);
P_r_UT(:, :, 1) = const.P0(1:2, 1:2);
P_v_UT(:, :, 1) = const.P0(3:4, 3:4);
P_r_UT(:, :, 1) = const.P0(1:2, 1:2);
P_v_UT(:, :, 1) = const.P0(3:4, 3:4);

eigen_UT_r = zeros(1, dt-1);
eigen_UT_v = zeros(1, dt-1);
eigen_UT_r(1) = 3*sqrt(max(eig(P_r_UT(:, :, 1))));
eigen_UT_v(1) = 3*sqrt(max(eig(P_v_UT(:, :, 1))));
 
for j = 2:dt

    for i = 1:2*n+1
        % Propagate the i-th sigma point using PCRTBP model
        [xf(:, i, j), ~, ~, xx, tt] = propagate_PCRTBP(const.ti, sigma_points(:, i), tspan(j), const.mu); 
    end

    mean_UT(:, j)=sum(W_m.*xf(:,:,j),2);
    P_UT(:,:,j)=W_c.*(xf(:,:,j)-mean_UT(:,j))*(xf(:,:,j)-mean_UT(:,j))';

     % Extract position and velocity covariance sub-matrices
    P_r_UT(:, :, j) = P_UT(1:2, 1:2, j);  % Position covariance
    P_v_UT(:, :, j) = P_UT(3:4, 3:4, j);  % Velocity covariance

    eigen_UT_r(j) = 3*sqrt(max(eig(P_r_UT(:, :, j))));
    eigen_UT_v(j) = 3*sqrt(max(eig(P_v_UT(:, :, j))));

end


for l = 5
    figure;
    hold on;
    plot(mean_LinCov(1, l), mean_LinCov(2, l), 'r.', 'MarkerSize', 10, 'LineWidth',2);  
    plot(mean_UT(1, l), mean_UT(2, l), 'b.', 'MarkerSize', 10, 'LineWidth',2);  
    error_ellipse(P_r_LinCov(:, :, l), mean_LinCov(1:2, l), 'conf', const.conf, 'style', 'r-', 'LineWidth', 1)
    error_ellipse(P_r_UT(:, :, l), mean_UT(1:2, l), 'conf', const.conf, 'style', 'b-', 'LineWidth', 1)
    xlabel('X');
    ylabel('Y');
    title(sprintf('Position LinCov vs UT @ t = %.4f', tspan(l)));  % Format time for title
    legend('Mean LinCov', 'Mean UT', 'LinCov Error Ellipse', 'UT Error Ellipse', 'Location', 'best')
    grid on
    hold off
end

for l = 5
    figure;
    hold on;
    plot(mean_LinCov(3, l), mean_LinCov(4, l), 'r.', 'MarkerSize', 10, 'LineWidth',2);  
    plot(mean_UT(3, l), mean_UT(4, l), 'b.', 'MarkerSize', 10, 'LineWidth',2);  
    error_ellipse(P_v_LinCov(:, :, l), mean_LinCov(3:4, l), 'conf', const.conf, 'style', 'r-', 'LineWidth', 1)
    error_ellipse(P_v_UT(:, :, l), mean_UT(3:4, l), 'conf', const.conf, 'style', 'b-', 'LineWidth', 1)
    xlabel('$\mathbf{V_x}$', 'Interpreter','latex');
    ylabel('$\mathbf{V_y}$', 'Interpreter','latex');
    title(sprintf('Velocity LinCov vs UT @ t = %.4f', tspan(l)));  % Format time for title
    legend('Mean LinCov', 'Mean UT', 'LinCov Error Ellipse', 'UT Error Ellipse', 'Location', 'best')
    grid on
    hold off
end


%% 2) Perform the same uncertainty propagation process on the same time grid using a Monte Carlo
   
clc; clearvars -except const mean_LinCov P_r_LinCov mean_UT P_r_UT eigen_LinCov_r eigen_LinCov_v eigen_UT_r eigen_UT_v  P_v_LinCov P_v_UT

dt = 5;
tspan = linspace(const.ti, const.tf, dt);
samples = 1000;

% Preallocate arrays
mean_MC = zeros(4, dt);    % Sample mean at each time step
P_MC = zeros(4, 4, dt);  % Sample covariance at each time step

P_r_MC = zeros(2, 2, dt);
P_v_MC = zeros(2, 2, dt);
P_MC(:, :, 1) = const.P0;
mean_MC(:, 1) = const.X0;
P_r_MC(:, :, 1) = const.P0(1:2, 1:2);
P_v_MC(:, :, 1) = const.P0(3:4, 3:4);
P_r_MC(:, :, 1) = const.P0(1:2, 1:2);
P_v_MC(:, :, 1) = const.P0(3:4, 3:4);

eigen_MC_r = zeros(1, dt);
eigen_MC_v = zeros(1, dt);
eigen_MC_r(1) = 3*sqrt(max(eig(P_r_MC(:, :, 1))));
eigen_MC_v(1) = 3*sqrt(max(eig(P_v_MC(:, :, 1))));

% Generate initial samples from the multivariate normal distribution
samples_MC = mvnrnd(const.X0, const.P0, samples)';

for j = 2:5
    for i = 1:samples
        [xf(:, i, j), ~, ~, xx, tt] = propagate_PCRTBP(const.ti, samples_MC(:, i), tspan(j), const.mu); 
    end
    
     % Compute sample mean
    mean_MC(:, j) = mean(xf(:, :, j), 2);

    % Compute sample covariance
    P_MC(:, :, j) = cov(xf(:, :, j)');  

     % Extract position and velocity covariance sub-matrices
    P_r_MC(:, :, j) = P_MC(1:2, 1:2, j);  % Position covariance
    P_v_MC(:, :, j) = P_MC(3:4, 3:4, j);  % Velocity covariance

    eigen_MC_r(j) = 3*sqrt(max(eig(P_r_MC(:, :, j))));
    eigen_MC_v(j) = 3*sqrt(max(eig(P_v_MC(:, :, j))));
    
end 

for l = 5
    figure;
    hold on;
    plot(xf(1, :, l), xf(2, :, l), '.', 'MarkerSize', 5, 'DisplayName', 'LinCov', 'LineWidth',2);  % Black dots for MC samples
    error_ellipse(P_r_LinCov(:, :, l), mean_LinCov(1:2, l), 'conf', const.conf, 'style', 'r-', 'LineWidth', 1)
    error_ellipse(P_r_UT(:, :, l), mean_UT(1:2, l), 'conf', const.conf, 'style', 'b-', 'LineWidth', 1)
    error_ellipse(P_r_MC(:, :, l), mean_MC(1:2, l), 'conf', const.conf, 'style', 'g-', 'LineWidth', 1)
    xlabel('X');
    ylabel('Y');
    title(sprintf('Position @ t = %.4f', tspan(l)), 'FontSize', 14);  % Format time for title
    legend('MC Samples', 'LinCov Error Ellipse', 'UT Error Ellipse', 'MC Error Ellipse','Location', 'best', 'FontSize', 12)
    grid on
    hold off
end

for l = 5
    figure;
    hold on;
    plot(xf(3, :, l), xf(4, :, l), '.', 'MarkerSize', 5, 'DisplayName', 'LinCov', 'LineWidth',2);  % Black dots for MC samples
    error_ellipse(P_v_LinCov(:, :, l), mean_LinCov(3:4, l), 'conf', const.conf, 'style', 'r-', 'LineWidth', 1)
    error_ellipse(P_v_UT(:, :, l), mean_UT(3:4, l), 'conf', const.conf, 'style', 'b-', 'LineWidth', 1)
    error_ellipse(P_v_MC(:, :, l), mean_MC(3:4, l), 'conf', const.conf, 'style', 'g-', 'LineWidth', 1)
    xlabel('VX');
    ylabel('VY');
    title(sprintf('Velocity @ t = %.4f', tspan(l)),'FontSize', 14);  % Format time for title
    legend('MC Samples', 'LinCov Error Ellipse', 'UT Error Ellipse', 'MC Error Ellipse','Location', 'best', 'FontSize', 12)
    grid on   
    hold off
end

figure
hold on;
% Plot the data with markers at each time point in tspan
plot(tspan, eigen_LinCov_r, '.', 'MarkerSize', 15, 'DisplayName', 'LinCov', 'LineWidth',2); 
plot(tspan, eigen_UT_r, '.', 'MarkerSize', 15, 'DisplayName', 'UT', 'LineWidth',2);  
plot(tspan, eigen_MC_r, '.', 'MarkerSize', 15, 'DisplayName', 'MC', 'LineWidth',2);  
% Display the legend
legend('Location', 'best', 'FontSize', 14)
% Enable grid for better visualization
grid on;
% Label the axes and set the title
xlabel('Time', 'Interpreter', 'none', 'FontSize', 14); % Default font
ylabel('Max Eigenvalue of Position Covariance', 'Interpreter', 'none', 'FontSize', 14); % Math symbols in LaTeX
title('Time Evolution of Max Eigenvalues of Position Covariance','FontSize', 14);
hold off

figure
hold on;
% Plot the data with markers at each time point in tspan
plot(tspan, eigen_LinCov_v, '.', 'MarkerSize', 15, 'DisplayName', 'LinCov', 'LineWidth',2); 
plot(tspan, eigen_UT_v, '.', 'MarkerSize', 15, 'DisplayName', 'UT', 'LineWidth',2);  
plot(tspan, eigen_MC_v, '.', 'MarkerSize', 15, 'DisplayName', 'MC', 'LineWidth',2);  
% Display the legend
legend('Location', 'best', 'FontSize', 14)
% Enable grid for better visualization
grid on;
% Label the axes and set the title
xlabel('Time', 'FontSize', 14);
ylabel('Max Eigenvalue of Velocity Covariance', 'FontSize', 14);
title('Time Evolution of Max Eigenvalues of Velocity Covariance', 'FontSize', 14);
hold off

% Loop over each state component

clc
ll = 5;
components = {'X', 'Y', 'V_{x}', 'V_{y}'};
for i = 1:4
    %subplot(2, 2, i);  % Arrange subplots
    %hold on
    figure
    qqplot(xf(i, :, ll));  % Generate Q-Q plot for component i
    title(['Q-Q plot for ', components{i}], 'FontSize', 14);   
    xlabel(['Theoretical ', components{i}], 'FontSize', 14);
    ylabel(['Real ', components{i}], 'FontSize', 14);
    grid on 
end

for i = 1:4
    [h(i), pValue(i)] = kstest(xf(i, :, ll), 'CDF', makedist('Normal', 'mu', mean(xf(i, :, ll)), 'sigma', std(xf(i, :, ll))));
        if h(i) == 0    
            disp('The sample is likely normally distributed.');
            disp(pValue(i));
        else
            disp('The sample is likely not normally distributed.');
            disp(pValue(i))
        end
    figure
    histogram(xf(i, :, ll), 20, 'Normalization', 'pdf');
    hold on;
    x = linspace(min(xf(i, :, ll)), max(xf(i, :, ll)), 100);
    y = normpdf(x, mean(xf(i, :, ll)), std(xf(i, :, ll)));
    plot(x, y, 'r', 'LineWidth', 2);
    title(sprintf('Histogram vs Normal Distribution (Component %d)', i));
    xlabel('State Component Value');
    ylabel('Probability Density');
end

% Save Figures

% saveAllFiguresToPath('/Users/pierazhi/Documents/MATLAB/GNC/Assignment 1 - Ex2/Images')


%% Functions

function h=error_ellipse(varargin)
% ERROR_ELLIPSE - plot an error ellipse, or ellipsoid, defining confidence region
%    ERROR_ELLIPSE(C22) - Given a 2x2 covariance matrix, plot the
%    associated error ellipse, at the origin. It returns a graphics handle
%    of the ellipse that was drawn.
%
%    ERROR_ELLIPSE(C33) - Given a 3x3 covariance matrix, plot the
%    associated error ellipsoid, at the origin, as well as its projections
%    onto the three axes. Returns a vector of 4 graphics handles, for the
%    three ellipses (in the X-Y, Y-Z, and Z-X planes, respectively) and for
%    the ellipsoid.
%
%    ERROR_ELLIPSE(C,MU) - Plot the ellipse, or ellipsoid, centered at MU,
%    a vector whose length should match that of C (which is 2x2 or 3x3).
%
%    ERROR_ELLIPSE(...,'Property1',Value1,'Name2',Value2,...) sets the
%    values of specified properties, including:
%      'C' - Alternate method of specifying the covariance matrix
%      'mu' - Alternate method of specifying the ellipse (-oid) center
%      'conf' - A value betwen 0 and 1 specifying the confidence interval.
%        the default is 0.5 which is the 50% error ellipse.
%      'scale' - Allow the plot the be scaled to difference units.
%      'style' - A plotting style used to format ellipses.
%      'clip' - specifies a clipping radius. Portions of the ellipse, -oid,
%        outside the radius will not be shown.
%
%    NOTES: C must be positive definite for this function to work properly.
default_properties = struct(...
  'C', [], ... % The covaraince matrix (required)
  'mu', [], ... % Center of ellipse (optional)
  'conf', 0.5, ... % Percent confidence/100
  'scale', 1, ... % Scale factor, e.g. 1e-3 to plot m as km
  'style', '', ...  % Plot style
  'clip', inf,...
  'LineWidth', ''); % Clipping radius
if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.C = varargin{1};
  varargin(1) = [];
end
if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.mu = varargin{1};
  varargin(1) = [];
end
if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.conf = varargin{1};
  varargin(1) = [];
end
if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.scale = varargin{1};
  varargin(1) = [];
end
if length(varargin) >= 1 & ~ischar(varargin{1})
  error('Invalid parameter/value pair arguments.') 
end
prop = getopt(default_properties, varargin{:});
C = prop.C;
if isempty(prop.mu)
  mu = zeros(length(C),1);
else
  mu = prop.mu;
end
conf = prop.conf;
scale = prop.scale;
style = prop.style;
if conf <= 0 | conf >= 1
  error('conf parameter must be in range 0 to 1, exclusive')
end
[r,c] = size(C);
if r ~= c | (r ~= 2 & r ~= 3)
  error(['Don''t know what to do with ',num2str(r),'x',num2str(c),' matrix'])
end
x0=mu(1);
y0=mu(2);
% Compute quantile for the desired percentile
k = sqrt(qchisq(conf,r)); % r is the number of dimensions (degrees of freedom)
hold_state = get(gca,'nextplot');
if r==3 & c==3
  z0=mu(3);
  
  % Make the matrix has positive eigenvalues - else it's not a valid covariance matrix!
  if any(eig(C) <=0)
    error('The covariance matrix must be positive definite (it has non-positive eigenvalues)')
  end
  % C is 3x3; extract the 2x2 matricies, and plot the associated error
  % ellipses. They are drawn in space, around the ellipsoid; it may be
  % preferable to draw them on the axes.
  Cxy = C(1:2,1:2);
  Cyz = C(2:3,2:3);
  Czx = C([3 1],[3 1]);
  [x,y,z] = getpoints(Cxy,prop.clip);
  h1=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
  [y,z,x] = getpoints(Cyz,prop.clip);
  h2=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
  [z,x,y] = getpoints(Czx,prop.clip);
  h3=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
  
  [eigvec,eigval] = eig(C);
  [X,Y,Z] = ellipsoid(0,0,0,1,1,1);
  XYZ = [X(:),Y(:),Z(:)]*sqrt(eigval)*eigvec';
  
  X(:) = scale*(k*XYZ(:,1)+x0);
  Y(:) = scale*(k*XYZ(:,2)+y0);
  Z(:) = scale*(k*XYZ(:,3)+z0);
  h4=surf(X,Y,Z);
  colormap gray
  alpha(0.3)
  camlight
  if nargout
    h=[h1 h2 h3 h4];
  end
elseif r==2 & c==2
  % Make the matrix has positive eigenvalues - else it's not a valid covariance matrix!
  if any(eig(C) <=0)
    error('The covariance matrix must be positive definite (it has non-positive eigenvalues)')
  end
  [x,y,z] = getpoints(C,prop.clip);
   h1 = plot(scale*(x0 + k*x), scale*(y0 + k*y), prop.style);
   set(h1, 'LineWidth', prop.LineWidth);  % Set LineWidth for the plot  set(h1,'zdata',z+1)
  if nargout
    h=h1;
  end
else
  error('C (covaraince matrix) must be specified as a 2x2 or 3x3 matrix)')
end
%axis equal
set(gca,'nextplot',hold_state);
%---------------------------------------------------------------
% getpoints - Generate x and y points that define an ellipse, given a 2x2
%   covariance matrix, C. z, if requested, is all zeros with same shape as
%   x and y.

end
%---------------------------------------------------------------
function [x,y,z] = getpoints(C,clipping_radius)
n=100; % Number of points around ellipse
p=0:pi/n:2*pi; % angles around a circle
[eigvec,eigval] = eig(C); % Compute eigen-stuff
xy = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
x = xy(:,1);
y = xy(:,2);
z = zeros(size(x));
% Clip data to a bounding radius
if nargin >= 2
  r = sqrt(sum(xy.^2,2)); % Euclidian distance (distance from center)
  x(r > clipping_radius) = nan;
  y(r > clipping_radius) = nan;
  z(r > clipping_radius) = nan;
end
end
%---------------------------------------------------------------
function x=qchisq(P,n)
% QCHISQ(P,N) - quantile of the chi-square distribution.
if nargin<2
  n=1;
end
s0 = P==0;
s1 = P==1;
s = P>0 & P<1;
x = 0.5*ones(size(P));
x(s0) = -inf;
x(s1) = inf;
x(~(s0|s1|s))=nan;
for ii=1:14
  dx = -(pchisq(x(s),n)-P(s))./dchisq(x(s),n);
  x(s) = x(s)+dx;
  if all(abs(dx) < 1e-6)
    break;
  end
end
end
%---------------------------------------------------------------
function F=pchisq(x,n)
% PCHISQ(X,N) - Probability function of the chi-square distribution.
if nargin<2
  n=1;
end
F=zeros(size(x));
if rem(n,2) == 0
  s = x>0;
  k = 0;
  for jj = 0:n/2-1;
    k = k + (x(s)/2).^jj/factorial(jj);
  end
  F(s) = 1-exp(-x(s)/2).*k;
else
  for ii=1:numel(x)
    if x(ii) > 0
      F(ii) = quadl(@dchisq,0,x(ii),1e-6,0,n);
    else
      F(ii) = 0;
    end
  end
end
%---------------------------------------------------------------
end
%---------------------------------------------------------------
function f=dchisq(x,n)
% DCHISQ(X,N) - Density function of the chi-square distribution.
if nargin<2
  n=1;
end
f=zeros(size(x));
s = x>=0;
f(s) = x(s).^(n/2-1).*exp(-x(s)/2)./(2^(n/2)*gamma(n/2));
end
%---------------------------------------------------------------
function properties = getopt(properties,varargin)
%GETOPT - Process paired optional arguments as 'prop1',val1,'prop2',val2,...
%
%   getopt(properties,varargin) returns a modified properties structure,
%   given an initial properties structure, and a list of paired arguments.
%   Each argumnet pair should be of the form property_name,val where
%   property_name is the name of one of the field in properties, and val is
%   the value to be assigned to that structure field.
%
%   No validation of the values is performed.
%
% EXAMPLE:
%   properties = struct('zoom',1.0,'aspect',1.0,'gamma',1.0,'file',[],'bg',[]);
%   properties = getopt(properties,'aspect',0.76,'file','mydata.dat')
% would return:
%   properties = 
%         zoom: 1
%       aspect: 0.7600
%        gamma: 1
%         file: 'mydata.dat'
%           bg: []
%
% Typical usage in a function:
%   properties = getopt(properties,varargin{:})
% Process the properties (optional input arguments)
prop_names = fieldnames(properties);
TargetField = [];
for ii=1:length(varargin)
  arg = varargin{ii};
  if isempty(TargetField)
    if ~ischar(arg)
      error('Propery names must be character strings');
    end
    f = find(strcmp(prop_names, arg));
    if length(f) == 0
      error('%s ',['invalid property ''',arg,'''; must be one of:'],prop_names{:});
    end
    TargetField = arg;
  else
    % properties.(TargetField) = arg; % Ver 6.5 and later only
    properties = setfield(properties, TargetField, arg); % Ver 6.1 friendly
    TargetField = '';
  end
end
if ~isempty(TargetField)
  error('Property names and values must be specified in pairs.');
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
