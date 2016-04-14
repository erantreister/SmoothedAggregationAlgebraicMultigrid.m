function [Asparsified_T] = SparsifyCollapsing(AT,PT,RT,P0T,R0T)

% tic
A0T = P0T*AT*R0T;
RP0T = P0T*RT;
R0P = (PT*R0T)';
A0 = A0T';

% SparsifyCollapsedGalerkinMex_new(A0T,RT,AT,PT,RP0T,R0P,A0); % this

% program includes multyplying the Galerkin inside.


% disp('Sparsening(Ac_omega,Ac_0)')
% disp('multiplying RAP')

SparsifyCollapsingMex(A0T,PT*AT*RT,RP0T,R0P,A0);


Asparsified_T = A0T;


if debuging()
    x = ones(size(A0T,2),1);
    if norm(A0' - A0,1)<1e-14 && norm(Asparsified_T' - Asparsified_T,1)>1e-10
        disp(['distance from symmetric Ac: ',num2str(norm(Asparsified_T' - Asparsified_T,1))]);
    end
    if (norm(Asparsified_T'*x - RT'*(AT'*(PT'*x))) > 1e-10)
        error('right null space not preserved');
    end
    if (norm(Asparsified_T*x - PT*(AT*(RT*x))) > 1e-10)
        error('left null space not preserved');
    end
end
%Asparsified_T = compareAndCorrect(Asparsified_T,AcT2);

return;
% figure;
% plot(eig(full(A0')),'*r');
% hold on;
% plot(eig(full(ActCOPY)),'ob');
% hold on;
% plot(eig(full(Asparsified_T )),'*g');
% hold on;
% legend('AGG','SA','SpSA');