function [x,y] = generateProtoypes(AT,x0,y0)
% clc
% AT = generateMatrix('ConvDiff2D',[256,256])';
% [x,S,y] = svds(AT,1,0);

n = size(AT,1);
x = CGNR(AT',x0,zeros(n,1),50);
y = CGNR(AT,y0,zeros(n,1),50);
x = x*(length(x)/norm(x));
y = y*(length(y)/norm(y));


% norm(AT'*x - y*S) /  norm(AT,1);
% S

% x = ones(size(AT,1),1);
% y = x;
% s = 1; %1./norm(AT,1);

% eig(full([zeros(size(AT)) , AT' ; AT , zeros(size(AT))]))
% for k = 1:50
% %     disp('-------------')
%     sigma = 0;
%     y = y - s*(AT*y  - sigma*x);
%     x = x - s*(AT'*x - sigma*y);
%     x = x/norm(x);
%     y = y/norm(y);
%     sigma = min(norm(AT'*x),norm(AT*y));
%     norm(AT'*x - sigma*y) /  norm(AT,1)
% end


function x = CGNR(A,x,b,nu) 
r = A'*(b - A*x);
w = r;
z = A'*(A*w);
a = (r'*w)/(w'*z);
x = x + a*w;
for i = 1:nu
    r = r - a*z;
    B = (r'*z)/(w'*z);
    w = r + B*w;
    z = A'*(A*w);
    a = (r'*w)/(w'*z);
    x = x + a*w;
end
return


return;
