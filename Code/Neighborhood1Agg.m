function Q = Neighborhood1Agg(S,order)
% This function builds an aggregation matrix Q according to the distance 
% one strength-based neighborhood aggregation algorithm. The input matrix S 
% is the symmetrized strength of connection matrix with 1's on its
% diagonal.

n = size(S,1);
% used = zeros(n,1);
% rowQ = zeros(n,1);
% colQ = zeros(n,1);
% 
% RCV = getRCV(S,0,0);
% % compute the partial aggregates Qhat_j
% numAggs = 0;
% prev = 1;
% next = 0;
% for i = 1:n
%     neighbors = RCV.Idx(RCV.Starts(i):RCV.Ends(i));
%     if sum(used(neighbors))==0
%         next = next+length(neighbors);
%         numAggs = numAggs+1;
%         rowQ(prev:next) = neighbors;
%         colQ(prev:next) = numAggs;
%         prev = next+1;
%         used(neighbors) = 1;
%     end
% end

% mean_card = nzmax(S)/size(S,2);
% card = sum(S~=0,2);
% card = card > 40*mean_card & card > 0.1*size(S,2);
% S(:,card) = 0;



if nargin == 1
    [numAggs,next,used,rowQ,colQ] = NeighborhoodAgg(S);
else
    [numAggs,next,used,rowQ,colQ] = NeighborhoodAggOrdered(S,order);
end


% the jth row of Qhat is the partial aggregate Qhat_j
Qhat = sparse(colQ(1:next),rowQ(1:next),ones(next,1),numAggs,n);

if next < n
    remain = find(used==0); % remaining unused points
    % columns of Nhat are strong neighborhoods of remaining unused points
    % and set all nonzero entries of Nhat equal to 1
    % Nhat = S(:,remain)>0;
	Nhat = S(:,remain)>0;
    % the ij entry of Qhat*Nhat is equal to card(Qhat_i intersect Nhat_j)
    [~,maxcard] = max(Qhat*Nhat,[],1);
    % update index arrays to obtain the full aggregates and then build Q
    rowQ((next+1):n,1) = remain;
    colQ((next+1):n,1) = maxcard';
    Q = sparse(rowQ,colQ,ones(n,1),n,numAggs);
else
    Q = Qhat'; 
end
if debuging()
    if sum(sum(abs(Q),1)==0)>0||sum(sum(abs(Q),2)==0)
        error('zero row or column in P0!!!');
    end
end


