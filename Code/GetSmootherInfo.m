function [info_rel] = GetSmootherInfo(params,AT)
info_rel = [];
info_rel.method = params.relax_method;
switch params.relax_method
    case {'L1'}
%         D = diag(AT);
%         omega = 1;
%         info_rel.Q = full(1./((1-omega)*D + omega*sum(abs(AT),1)'));
%         info_rel.Q = full(1./(sum(abs(AT))'));
        info_rel.Q = GetDiagPreconditioner(AT,'1');
    case {'SPAI','GSSPAI'}
%         D = diag(AT);
%         info_rel.Q = full(D./(sum(AT.^2,1)'));
        info_rel.Q = GetDiagPreconditioner(AT,'S');
    case {'JAC','GS'} 
        info_rel.Q = full(1./diag(AT));
    case 'MultiColor'
        info_rel.Q = full(1./diag(AT));
        [info_rel.multiColorC,info_rel.multiColorI] = multiColorSetup(AT,[]);
        I = info_rel.multiColorI;
        disp(['number of colors: ',num2str(length(I)-1),' Distributed: ' num2str(I(2:end)-I(1:end-1)) ]);
    otherwise
        error('Unsupported P smoothing');
end
% Determine parameters for relaxations
if strcmp(params.RelParamType,'fixed')
    info_rel.omega = params.RelParam;
    if strcmp(info_rel.method,'GS')==1 || strcmp(info_rel.method,'GSSPAI')==1 
        info_rel.Q = info_rel.Q.*params.RelParam;
    end
elseif strcmp(params.RelParamType,'minFRO')
    Qrel = info_rel.Q;
    info_rel.omega = (diag(AT)'*Qrel)/(sum(((Qrel.^2)').*sum(AT.^2,1)));
else
    error('omega_relaxation - not determined')
end

