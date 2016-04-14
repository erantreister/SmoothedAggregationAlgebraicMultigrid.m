function [x,r] = Jacobi(A,x,r,omega,nu,b)
Q = 1./diag(A);
if nargout == 2
    [x,r] = Relax(A,x,r,omega,nu,Q,b);
else
    x = Relax(A,x,r,omega,nu,Q,b);
end

return;