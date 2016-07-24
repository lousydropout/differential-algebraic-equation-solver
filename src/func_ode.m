function result = func_ode(w, op, param, ic)  % function for ODE
  a = param.g / param.l;
  result = op.D2*w + a*sin(w);
  result(1) = w(1) - ic.theta0;
  result(2) = op.D(1,:) * w - ic.omega0;
endfunction
