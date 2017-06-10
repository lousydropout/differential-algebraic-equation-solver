clc; clear; close;

N = 50;     % no. collocation points
t_start = 0; 
t_finish = 10;
t_step = 3.1;
[op.D, op.t] = cheb(N-1, 0, t_step);  % set up grid & derivative op.


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
t     = t_start;


%% Create vectors to store solutions
sol.t = [];
sol.x = [];
sol.u = [];
sol.y = [];
sol.v = [];
sol.f = [];


while (t < t_finish)
    
    if (t + t_step > t_finish)
        tmp    = t_finish - t;
        op.D   = op.D * t_step / tmp;
        op.t   = op.t / t_step * tmp;
        t_step = tmp;
    endif
    
    %% set up vectors
    % used for DAE case
    x = ic.x0 * ones(N, 1);
    u = ic.u0 * ones(N, 1);
    y = ic.y0 * ones(N, 1);
    v = ic.v0 * ones(N, 1);
    f = ic.f0 * ones(N, 1);
    w = [x; u; y; v; f];


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

    
    %% Store solutions
    sol.t = [sol.t; t + op.t];
    sol.y = [sol.y; y];
    

    %% Update IC for new time step
    ic.x0 = x (N);
    ic.u0 = u (N);
    ic.y0 = y (N);
    ic.v0 = v (N);
    ic.f0 = f (N);
    t     = t + t_step;
    
endwhile    

%% Plot and compare y for ODE and DAE solutions
hold off; plot(sol.t, sol.y); 