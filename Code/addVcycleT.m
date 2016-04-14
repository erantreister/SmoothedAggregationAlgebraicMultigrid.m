function x = addVcycleT(x,r,b,params,MGOP,curr_level)
% if (sum(isnan(x))>0)||(sum(isnan(r))>0)||(sum(isnan(b))>0)
%     error('ERAN: x is nan');
% end
max_levels = length(MGOP)-1;
AT          = MGOP{curr_level}.AT;
R          = MGOP{curr_level}.R;
P          = MGOP{curr_level}.P;
info_rel       = MGOP{curr_level}.info_rel;

if params.nu1==1 || (curr_level==1 && params.oneRelaxOnFineGrid)
    [x,r] = RelaxT(AT,x,b,r,1,info_rel,'F');
elseif params.nu1==2
    [x,r] = RelaxT(AT,x,b,r,1,info_rel,'S');
end

% if (sum(isnan(x))>0)
%     error('ERAN: x is nan');
% end



% [x,r] = RelaxT(AT,x,b,r,params.nu1,info_rel,'B');
b_c = R*r;
if curr_level == max_levels
    e_c = solveCoarsest(MGOP{curr_level+1}.AT',b_c);
else
    r_c = -b_c;
    x_c = zeros(size(r_c));
    e_c = addVcycleT(x_c,r_c,b_c,params,MGOP,curr_level+1);
    if strcmp(params.cycle_type,'W')
        e_c = addVcycleT(e_c,[],b_c,params,MGOP,curr_level+1);
    end
    if strcmp(params.cycle_type,'F')
        params.cycle_type = 'V';
        e_c = addVcycleT(e_c,[],b_c,params,MGOP,curr_level+1);
    end
end
e = P*e_c;
    
if params.nu2==1 || (curr_level==1 && params.oneRelaxOnFineGrid)
    e = RelaxT(AT,e,r,[],1,info_rel,'B');
elseif params.nu2==2
    e = RelaxT(AT,e,r,[],1,info_rel,'S');
end

if strcmp(params.OC,'rec_minres')==1
    if (curr_level > 1)
        Ae = AT'*e;
%         alpha = (r'*Ae)/(Ae'*Ae);
        alpha = (r'*e)/(e'*Ae);
        %     r = r - alpha*Ae;
        x = x - alpha*e;
    else
        x = x - e;
    end
else
    if  params.OC~=1
        e = params.OC*e;
    end
    x = x - e;
end
return;

function [x] = solveCoarsest(A,b)
% [U,S,V] = svd(full(A));
% Sinv = 1./(diag(S)+eps);
% if Sinv(end)/norm(A,1) > 1e+10
%     Sinv(end) = 0;
% end
% % disp('the null-free solution')
% x = V*(Sinv.*(U'*b));

% x = A\b;

x = (A + 1e-14*norm(A,1)*speye(size(A)))\b;
removeConstant = sum(abs(A*ones(size(b))))/size(A,2) < 1e-11;
if removeConstant
%     disp('removing const');
    x = x - sum(x)/length(x); % remove constant from solution
end

if (sum(isnan(x))>0)
    error('ERAN: x is nan');
end
return;
