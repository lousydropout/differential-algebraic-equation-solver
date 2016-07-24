function Jac = jacobian_dae(w, op, param)  % Jacobian for DAE
  N = param.N;
  x = w(      1:   N); X = diag(x);
  y = w(2*N + 1: 3*N); Y = diag(y);
  f = w(4*N + 1: 5*N); F = diag(f);

  l  = param.l; ml = param.m * l; mlg = ml * param.g;
  D = op.D;

  m1 = -eye(N); zero = zeros(N, N); 
  
  Jac = [D,    m1,   zero, zero, zero;
  	 F,    ml*D, zero, zero, X;
  	 zero, zero, D,    m1,   zero;
  	 zero, zero, F,    ml*D, Y;
  	 X,    zero, Y,    zero, zero];
  Jac(      1, :) = [1, zeros(1, 5*N - 1)];
  Jac(  N + 1, :) = [zeros(1,   N), 1, zeros(1, 4*N - 1)];
  Jac(2*N + 1, :) = [zeros(1, 2*N), 1, zeros(1, 3*N - 1)];
  Jac(3*N + 1, :) = [zeros(1, 3*N), 1, zeros(1, 2*N - 1)];
  Jac(4*N + 1, :) = [zeros(1, 4*N), 1, zeros(1,   N - 1)];

  Jac = sparse(Jac);
endfunction
