function [x,rnorms,r] = PCG(AT,x,r,max_iter,M,tol)
z = M(r);
j = 1;
rnorms = zeros(1,max_iter);
rnorms(1) = sqrt(r'*r);
ztr = z'*r;
p = z;
while (1)
    z = AT'*p; % here z is actually Ap
    alpha = ztr/(p'*z);
    x = x + alpha*p;
    r = r - alpha*z;
    rnorms(j) = sqrt(r'*r);
    if (rnorms(j) < tol) || (j >= max_iter)
        rnorms = rnorms(1:j);
        break;
    end
    z = M(r);
    ztr_new = (z'*r);
    beta = ztr_new/ztr;
    p = z + beta*p;
    ztr = ztr_new;
    j = j + 1;
end
return

% function [x,rnorms,rnew] = PCG(AT,x0,r0,max_iter,M,tol)
% z = M(r0);
% r = r0;
% p = z;
% j = 1;
% x = zeros(size(x0));
% rnorms = zeros(1,max_iter);
% rnorms(1) = sqrt(r0'*r0);
% while (1)
%     Ap = AT'*p;
%     alpha = (z'*r)/(p'*Ap);
%     x = x + alpha*p;
%     rnew = r - alpha*Ap;
%     rnorms(j) = sqrt(rnew'*rnew);
%     if (rnorms(j) < tol) || (j >= max_iter)
%         rnorms = rnorms(1:j);
%         break;
%     end
%     znew = M(rnew);
%     beta = (znew'*rnew)/(z'*r);
%     p = znew + beta*p;
%     r = rnew;
%     z = znew;
%     j = j + 1;
% end
% x = x0 + x;
% return
