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
% % For more information see papers:
% (1)Eran Treister and Irad Yavneh.
%    Square and Stretch Multigrid for Stochastic Matrix Eigenproblems.
%    Numerical Linear Algebra with Application, 17:229-251, 2010
% 
% (2)Eran Treister and Irad Yavneh.
%    On-the-fly Adaptive Smoothed Aggregation Multigrid for Markov Chains.
%    SIAM Journal on Scientific Computing (SISC), 33:2927-2949, 2011

function [aggr,M] = getAggregation(AT,x,params)
% A - sparese M-matrix
% x - positive right null-space prototype vector.
% method - currently only Bottom-Up is implemented. method should be
% ['BUp',num2str(size)] where size is an average support size integer.
% aggr - a map from each grid point 1..n to its aggregate's root node.
%        The indices that appead in the aggr array are the root nodes, and they represent 
%        the coarse level variables.
if (sum(isnan(x))>0)
    error('x is nan');
end
method = params.coarsening_method;

M=[];
if strcmp(method(1:end-1),'BUp')
    x = abs(x);
    thetta = 0.25;
    AggrSize = round(str2double(method(end)));
    n = length(x);
    [C,R,V] = find(AT);
    V = V'; C = int32(C'); 
    R = ToRowStarts(R,n,length(R));
    [starts,C,V] = WeightMatrix(R,C,V,n,length(V),x,thetta);
    aggr = BUpAgg(AggrSize,starts,C,V,n,length(V));
    aggr = aggrArray2P(aggr);
end


if strcmp(method,'N1')    
    M = StrengthMatrix(AT);
    if debuging()
        if sum(sum(abs(M),1)==0)>0||sum(sum(abs(M),2)==0)
            error('zero row or column in M!!!');
        end
    end
    S = (M >= params.P_gammaStrong);
    S = 0.5*(S + S');
    if debuging()
        if sum(sum(abs(S),1)==0)>0||sum(sum(abs(S),2)==0)
            error('zero row or column in S!!!');
        end
    end
    aggr = Neighborhood1Agg(S);    
end
if strcmp(method,'NN1')    
    M = StrengthMatrix(AT);
    if debuging()
        if sum(sum(abs(M),1)==0)>0||sum(sum(abs(M),2)==0)
            error('zero row or column in M!!!');
        end
    end
    S = (M>=params.P_gammaStrong);
    S = 0.5*(S + S');
    if debuging()
        if sum(sum(abs(S),1)==0)>0||sum(sum(abs(S),2)==0)
            error('zero row or column in S!!!');
        end
    end
    aggr = NeighborhoodAggNew(S);
    aggr = aggrArray2P(aggr);
end


if strcmp(method,'ON1')
    M = StrengthMatrix(AT);
    S = (M>=params.P_gammaStrong);
    S = 0.5*(S + S');
    
    %     BT = abs(AT);
    %     BT = 0.5*(BT+BT');
    %     D = toSparseDiag(1./sqrt(diag(BT)));
    %     BT = D*BT*D;
    %     S = BT.*(BT > 0.02);
    numNei = sum(S);
    n = length(numNei);
    [tmp,BupOrder] = sort(numNei + (0:n-1)/n,'ascend');
    aggr = Neighborhood1Agg(S,BupOrder);
    %     size(aggr,2)
end

return;


function [P] = aggrArray2P(aggr)
fine2coarse = zeros(size(aggr));
n = length(aggr);
I = find(aggr == 1:n);
fine2coarse(I) = 1:length(I);
aggr = double(fine2coarse(aggr));
if sum(aggr==0)>0
    error('nodes without aggregates');
end
P = sparse(1:n,aggr,ones(1,n));
return;