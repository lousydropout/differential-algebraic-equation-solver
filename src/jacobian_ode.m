function result = jacobian_ode(w, op, param)   % Jacobian for ODE
  N = param.N;
  a = param.g / param.l;
  result = op.D2 + diag(a*cos(w));
  result = [1 zeros(1, N-1);
	    op.D(1,:);
	    result(3:N, :)];
endfunction
