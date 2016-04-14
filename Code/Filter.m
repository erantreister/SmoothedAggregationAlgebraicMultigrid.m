function [AT_filtered] = Filter(AT,epsilon,M)
%FILTER Summary of this function goes here
%   Detailed explanation goes here
% AT_filtered = AT;
% return;
% D = toSparseDiag(1./sqrt(diag(AT)));
% Ind = abs(D*AT*D) > epsilon;



if isempty(M)
   M = StrengthMatrix(AT);
   
end
Ind = abs(M) > epsilon;


% invD = full(1./diag(AT));
% S = abs(AT);
% DiagOperate(S,invD,'L')
% Ind = S > epsilon;

% mean_card = nzmax(M)/size(M,2);
% card = sum(spones(M));
% card = card > 40*mean_card & card > 0.1*size(M,2);
% Ind(:,card) = 0;



% n = size(AT,2);
% Ind = abs(AT) - spdiags(diag(AT),0,n,n);
% m = max(Ind,[],2);
% Ind = spdiags(1./m + 1e-14,0,n,n)*Ind;
% Ind = (Ind>=5*epsilon) + speye(n);


AT_filtered = AT.*Ind;

o = ones(size(AT,1),1);
DiagOperate(AT_filtered,AT'*o - AT_filtered'*o,'A');
% AT_filtered(:,card) = 0;
return;
