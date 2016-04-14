% ---------------------------------------------------------------------------
% Copyright (2012): Eran Treister and Irad Yavneh. 
% This code is distributed under the terms of the GNU General Public
% License 2.0.
% 
% Permission to use, copy, modify, and distribute this software
% for any purpose without fee is hereby granted, provided that
% this entire notice is included in all copies of any software
% which is or includes a copy or modification of this software
% and in all copies of the supporting documentation for such
% software. This software is being provided "as is", without any
% express or implied warranty. In particular, the authors do not
% make any representation or warranty of any kind concerning the
% merchantability of this software or its fitness for any
% particular purpose."
% ---------------------------------------------------------------------------
function [P,R,xc,yc] = getAggOperators(aggr,x,y)
% This function generates aggregation operators P and R with x and y in
% their range, respectively. This version uses an l-2 normalization for
% each aggregate, s.t. P'*P = I and R*R' = I.
n = length(x);
P_bar = aggr;
% define adaptive P and R:
xc = sqrt(P_bar'*(x.*x));
P = toSparseDiag(x)*P_bar*toSparseDiag(1./xc);

% t = x./xc(I);
% P = sparse(1:n,I,t);
if isempty(y)
    R = P';
else
    yc = sqrt(P_bar'*(y.*y));
    R = toSparseDiag(1./yc)*P_bar'*toSparseDiag(y);
%     t = y./yc(I);
%     R = sparse(I,1:n,t);
end
% P = P>0;
% R = R>0;
return;
function Q = toSparseDiag(q)
n = length(q);
Q = spdiags(q,0,n,n);
return;