function [x,MGOP,trace_info] = SetupVcycleT(AT,x0,params,MGOP)
% nu = 0;

if isempty(MGOP)
    MGOP = cell(12,1);
end
MGOP{1}.AT = AT;
lev = 1;
n = size(AT,2);
NNZ = 0;
x = x0;
params.buildAggr = 1;
MGOP{1}.AT = AT;

while exactSolveInd(n,lev,params)==0
    AT = MGOP{lev}.AT;
    zero_vec = zeros(size(AT,2),1);
    info_rel = GetSmootherInfo(params,AT);
    NNZ = NNZ + nzmax(AT);
    
    if lev == 1 && isempty(x)
        x = ones(size(zero_vec));     
    end
    if params.isSymmetric==0
        if lev==1
            y = ones(size(zero_vec));
        end
    end
    
    if params.buildAggr
        x = ones(size(AT,1),1);
        [MGOP{lev}.aggr,M] = getAggregation(AT,x,params);
%         disp(['Coarsening: n->nc: ',num2str(size(MGOP{lev}.aggr,1)),'->',num2str(size(MGOP{lev}.aggr,2))])
    else
        M=[];
    end
    P0 = MGOP{lev}.aggr;    R0 = P0';
    P0 = P0'; R0 = R0';
    
    if strcmp(params.setup_method,'SA')||strcmp(params.setup_method,'hSA')||strcmp(params.setup_method,'SpSA')
        
        if params.P_epsilonFilter~=0
            AT_F = Filter(AT,params.P_epsilonFilter,M);
        else
            AT_F = AT;
        end
        [Qcg,omega_cg] = smoothProlongationInfo(params,AT_F,P0);
%         D = toSparseDiag(omega_cg*Qcg);
%         P = P0 - (P0*AT_F)*D;        
        T = P0*AT_F;
        DiagOperate(T,omega_cg*Qcg,'L');
        P = P0 - T;
        if debuging()
            if sum(sum(abs(P),1)==0)>0||sum(sum(abs(P),2)==0)
                error('zero row or column in P!!!');
            end
        end
        if params.isSymmetric==0
%             R = R0 - D*(AT_F*R0);
            
            T = AT_F*R0;
            DiagOperate(T,omega_cg*Qcg,'R');
            R = R0 - T;
            if debuging()
                if sum(sum(abs(R),1)==0)>0||sum(sum(abs(R),2)==0)
                    error('zero row or column in P!!!');
                end
            end
        else
            R = P';
        end
    else
        P = P0; R = R0;
    end
%     figure;
%     spy(AT);
%     title(['Level: ',num2str(lev), ', size: ',num2str(n),'x',num2str(n),' #nz: ',num2str(nnz(AT))],'Fontsize',12);
    
    if strcmp(params.setup_method,'SA')||strcmp(params.setup_method,'AGG')
        
        Act = P*(AT*R);
        
        if debuging()
            if sum(sum(abs(AT),1)==0)>0||sum(sum(abs(AT),2)==0)
                error('zero row!!!');
            end
        end
    elseif strcmp(params.setup_method,'SpSA')
        
        Act = SparsifyCollapsing(AT,P,R,P0,R0);
    elseif strcmp(params.setup_method,'hSA')   
        Act = P*(AT*R0);
    end

    MGOP{lev}.info_rel = info_rel;
    MGOP{lev}.R = R';
    MGOP{lev}.P = P';
    MGOP{lev+1}.AT = Act;
    n = size(R,2);
    lev = lev+1;
end
trace_info.NNZ = NNZ + nzmax(MGOP{lev}.AT);
trace_info.levels = lev;
MGOP = MGOP(1:lev);
% figure;
% spy(Act);
% title(['Level: ',num2str(lev), ', size: ',num2str(n),'x',num2str(n),' #nz: ',num2str(nnz(Act))],'Fontsize',12);

return;
