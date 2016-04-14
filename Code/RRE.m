function [s] = RRE(X)
k = size(X,2)-2;
if k == 0
    s = X(:,1);
    return
end
if k < 0 
    error('too few vectors to MPE');
end
U = X(:,1:end-1);% size (N,k+1)
for i=1:(k+1)
    U(:,i) = X(:,i+1) - X(:,i);
end
[Q_k R_k] = qr(U,0);    % [N k+1] [k+1 k+1] 
e = ones(k+1,1);
d = R_k\(R_k'\e);
gamma = d/sum(d);

R_km1 = R_k(1:k , 1:k);
xi = zeros(k,1);
xi(1) = 1-gamma(1);     
for i=2:length(xi)
    xi(i) = xi(i-1) - gamma(i);
end
s = X(:,1) + Q_k(:,1:end-1)*(R_km1*xi);
