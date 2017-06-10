%% This is my Newton-Raphson solver
function [u, success] = solve(func, jacobian, w, op, param, ic)
  res = 1; iter = 1;
  u   = w;
  F   = func(u, op, param, ic);
  
  while (res > param.tol)
    % check to see if solution is blowing up
    if ((iter > 100) && (res > 1))
       "*** Newton-Raphson subroutine FAILED to converge ***"
       success = false;
       return;
    endif
    % Newton-Raphson iteration
    Jac = jacobian(u, op, param);
    du  = Jac \ F;
    u   = u - du;
    F   = func(u, op, param, ic);
    iter = iter + 1;
    res = norm(F);
  endwhile
  success = true;
endfunction
