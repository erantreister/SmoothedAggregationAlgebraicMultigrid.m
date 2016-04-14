function [x,rnorms,r] = GMRES(AT,x0,r0,inner,prec,TOL)
% if inner==1
%     x = MR(AT,r0,x0,prec);
%     return;
% end

% r0 = b-AT'*x0;
betta = sqrt(r0'*r0);
H = zeros(inner+1,inner);

Z = zeros(length(r0),inner);
V = zeros(length(r0),inner);
AZ = zeros(length(r0),inner);

w = r0/betta;
 

xi = zeros(inner+1,1);
xi(1) = betta;

rnorms = zeros(1,inner);

for j = 1:inner
    V(:,j) = w;    
    z = prec(w);
    w = AT'*z;
    AZ(:,j) = w;
    % GS:
    betta = V'*w;
    H(1:inner,j) = betta;
    w = w - V*betta;
%     % MGS:
%     for i=1:j
%         H(i,j) = w'*V(:,i);
%         w = w - H(i,j)*V(:,i);
%     end
%     toc
    Z(:,j) = z;
    betta = sqrt(w'*w);
    H(j+1,j) = betta;
    w = w*(1/betta);

%     y = (H(1:j+1,1:j)\xi(1:j+1));
%     norm(H(1:j+1,1:j)*y - xi(1:j+1))
    [Q,~] = qr(H(1:j+1,1:j));
    rnorms(j) = abs(Q(:,end)'*xi(1:j+1));
    if rnorms(j) < TOL
        rnorms = rnorms(1:j);
        break;
%         y = zeros(inner,1);
%         y(1:j) = (H(1:j+1,1:j)\xi(1:j+1));
%         x = x0 + Z*y; % multiplying the zeros is cheaper...go figure..
%         r = r0 - AZ*y;
%         return;
    end
end
y = pinv(H)*xi;
% y = (H\xi);
% norm(H*y - xi)
x = x0 + Z*y;
r = r0 - AZ*y;
return;

% function [x] = MR(AT,r0,x0,prec)
% e0 = prec(r0);
% Ae0 = AT'*e0;
% % alpha = 1;
% alpha = (Ae0'*r0)/(Ae0'*Ae0);
% x = x0 + alpha*e0;
% return;

