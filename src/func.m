function result = func (w, op, param, ic) % function for DAE
  N = param.N;
  x = w (      1:   N);
  u = w (  N + 1: 2*N);
  y = w (2*N + 1: 3*N);
  v = w (3*N + 1: 4*N);
  f = w (4*N + 1: 5*N);

  l  = param.l;
  ml = param.m * l;
  mlg = ml * param.g;
  D = op.D;

  % set up function
  result = zeros (5 * N, 1);
  result (      1:   N) = D * x - u;
  result (  N + 1: 2*N) = ml * (D*u) +  (f.*x);
  result (2*N + 1: 3*N) = D * y - v;
  result (3*N + 1: 4*N) = ml * (D*v) + (f.*y) - mlg;
  result (4*N + 1: 5*N) = 0.5 * (x.^2 + y.^2 - l^2);

  % correct function
  result (      1) = x (1) - ic.x0;
  result (  N + 1) = u (1) - ic.u0;
  result (2*N + 1) = y (1) - ic.y0;
  result (3*N + 1) = v (1) - ic.v0;
  result (4*N + 1) = f (1) - ic.f0;
endfunction

