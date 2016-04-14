function demo()
addpath('Code');
addpath('Code/MEXfunc');

% FORMING 2D Laplacian operator:
A = delsq(numgrid('S',513));
b1 = rand(size(A,2),1);
AT = A';clear A;


[x1,setup_info,results] = SolveLinearSystem(AT,b1);
disp(['~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~']);
disp(['~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~']);
disp('Relative norm of the residual:')
norm(AT'*x1-b1)/norm(b1)

rmpath('Code');
rmpath('Code\MEXfunc');
return;