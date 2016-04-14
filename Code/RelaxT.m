function [x,r] = RelaxT(AT,x,b,r,nu,info_rel,direction)
if nu == 0 
    return;
end
if strcmp(info_rel.method,'GS')==1 || strcmp(info_rel.method,'GSSPAI')==1 
    for i=1:nu
        x = GaussSeidelSMatTranspose(AT,x,b,info_rel.Q,direction);
    end    
elseif strcmp(info_rel.method,'MultiColor')==1
    x = MulticolorGaussSeidelSMatTranspose(AT,x,b,info_rel.Q,nu,info_rel.multiColorC,info_rel.multiColorI);
else
    omega = info_rel.omega;
    Q = info_rel.Q;
    if isempty(r)
        r = AT'*x-b;
    end
    for i=1:nu-1
        x = x - omega*(Q.*r);
        r = AT'*x-b;
    end
    x = x - omega*(Q.*r);
end

if nargout==2
    r = AT'*x-b;
end
return;