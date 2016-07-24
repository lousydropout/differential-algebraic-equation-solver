clc; clear; close;

N = 100;     % no. collocation points
t_start = 0; t_finish = 30;
[op.D, t] = cheb(N-1, t_start, t_finish);  % set up grid & derivative op.
op.D2 = op.D^2;  % Laplacian operator


%% parameters, held in scalar data structure 'param'
param.N   = N;     % no. collocation points
param.m   = 1;     % mass
param.l   = 13.7;    % length
param.g   = 10;    % gravitational acceleration
param.tol = 1.e-6;  % error tolerance


%% initial conditions, held in scalar data structure 'ic'
% used for DAE case
ic.x0 = 11;
ic.y0 = sqrt(param.l^2 - ic.x0^2);   
ic.u0 = 0;
ic.v0 = 0;
ic.f0 = param.m * param.g * ic.y0 / param.l; % initial tension
% used for ODE case
ic.theta0 = atan(ic.x0 / ic.y0); % initial angle
ic.omega0 = 0; % initial angular velocity


%% set up vectors
% used for DAE case
x = ic.x0 * ones(N, 1);
u = ic.u0 * ones(N, 1);
y = ic.y0 * ones(N, 1);
v = ic.v0 * ones(N, 1);
f = ic.f0 * ones(N, 1);
w = [x; u; y; v; f];
% used for ODE case
th = ic.theta0 * ones(N, 1);


%% Solve DAE
[u, success] = solve(@func_dae, @jacobian_dae, w, op, param, ic);
if (success == true)
  w = u;
  x = w(1:N);
  u = w(  N + 1: 2*N);
  y = w(2*N + 1: 3*N);
  v = w(3*N + 1: 4*N);
  f = w(4*N + 1: 5*N);
endif


%% Solve ODE
[u, success] = solve(@func_ode, @jacobian_ode, th, op, param, ic);
if (success == true)
   th = u;
   x_ode = param.l * sin(th);
   y_ode = sqrt(param.l^2 - x_ode.^2);
endif;


%% Plot and compare y for ODE and DAE solutions
hold off;plot(t, y); hold on; plot(t, y_ode, 'ro'); hold off;