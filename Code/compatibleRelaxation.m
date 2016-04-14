function [C] = compatibleRelaxation(AT)
n = size(AT,2);
Os = zeros(n,1);

rel_params.relax_method = 'SPAI';
rel_params.RelParamType = 'fixed';
rel_params.RelParam = 1;
info_rel = GetSmootherInfo(rel_params,AT);
e_0 = rand(n,1);
e_0 = RelaxT(AT,e_0,Os,[],7,info_rel,'F');
e_0 = RelaxT(AT,e_0,Os,[],7,info_rel,'B');
e_0 = abs(e_0);
rel_params.relax_method = 'JAC';
rel_params.RelParamType = 'fixed';
rel_params.RelParam = 1;
info_rel = GetSmootherInfo(rel_params,AT);
nu_cr = 1;
theta_cr = 0.7;

C = zeros(n,1)>0;
[rho,ratios] = iterateCR(AT,e_0,Os,nu_cr,info_rel,C);

M2 = abs(AT);
M2 = (M2*M2)';
while rho > theta_cr
    [C] = getCoarsePoints(AT,M2,ratios,theta_cr,C);
%     disp(['|C|/n:',num2str(sum(C)/n)])
    for i = 2
        [rho,ratios] = iterateCR(AT,e_0,Os,nu_cr,info_rel,C);
    end
%     disp(['rho:',num2str(rho)])
%     disp('*****************************');
end
return;

function [C] = getCoarsePoints(AT,M2,ratios,theta_cr,C)
k = 1;
n = size(AT,1);
[s_ratios,i_ratios] = sort(ratios,'descend');
while k <= n
    idx = i_ratios(k);
    if ratios(idx) > theta_cr
        C(idx) = 1; % disp(['choosing: ',num2str(idx)]);
        ratios(M2(:,idx)~=0) = 0;
    elseif (s_ratios(k) < theta_cr)
        break;
    end
    k = k+1;
end
C = C>0;
return;

function [rho,ratios] = iterateCR(AT,e,Os,nu_cr,info_rel,C)
% e - positive smooth error vector
e_old = e;
e(C) = 0;
for k=1:nu_cr
    e = RelaxT(AT,e,Os,[],1,info_rel,'F');
    e(C) = 0;
    e = RelaxT(AT,e,Os,[],1,info_rel,'B');
    e(C) = 0;
end
e_new = e;
ratios = abs(e_new)./(e_old + 1e-14);
rho = norm(e_new) / norm(e_old(~C));
return;


% function [N_ind] = getNeighborhood2(A,i)
% % N_ind = (A*A(:,i))~=0;
% N_ind = A(:,i)~=0;
% N_ind(i) = 0;
% N = find(N_ind);
% for k = 1:length(N)
%     N_ind = N_ind | A(:,N(k))~=0;
% end
% return;


% function [N] = getNeighborhood1(A,i)
% N = A(:,i)~=0;
% return;
% 
% function [N] = getNeighborhood(A,i,dist)
% if dist == 1
%     N = find(A(:,i)~=0);
% else
%     N_ind = A(:,i)~=0;
%     N = find(N_ind);
%     for k = 1:length(N)
%         N_ind(getNeighborhood(A,N(k),dist-1)) = 1;
%     end
%     N = find(N_ind);
% end
% return;
