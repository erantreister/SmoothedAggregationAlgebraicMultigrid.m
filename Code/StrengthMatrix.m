function S = StrengthMatrix(AT)
% n = size(AT,2);
% D = spdiags(diag(AT),0,n,n);
% S = -AT + D;
% m_rows = max(S)';
% M1 = S*spdiags(1./m_rows,0,n,n);
% M = M1 + speye(n);


n = size(AT,2);
S = (-AT);
DiagOperate(S,zeros(n,1),'E');
m_rows = max(S)';
m_rows = m_rows + 1e-16*max(m_rows);
if debuging()
    if sum(m_rows==0)>1
        error('I dont know..');
    end
end
DiagOperate(S,1./m_rows,'L');
DiagOperate(S,ones(n,1),'A');

% norm(S-M)
end

