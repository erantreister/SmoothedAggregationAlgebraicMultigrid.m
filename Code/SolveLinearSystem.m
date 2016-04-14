function [x,setup_info,results] = SolveLinearSystem(AT,b,x0,isSymmetric,setup_info,maxit,TOL,params,tracefile)
% SolveLinearSystem 
% Solves a linear system using aggregation based multigrid.
% Input - see below. The only requirements are AT, b. Note: AT is the
% transpose of the matrix to solve, A*x = b.
% Output:
% x - the solution.
% setup_info - a stuct with the MG setup. Can be reused as input for
% solving the same system with more RHSs.
% results - information on the run.

addpath('./Code/MEXfunc');
%% INPUT:

% AT - the transpose of the linear system's matrix.
% b - rhs.

if nargin < 4 || isempty(isSymmetric)
    if max(max(abs(AT-AT')))<1e-15
        disp('Symmetric');
        isSymmetric = 1;
    else
        disp('Non-Symmetric');
        isSymmetric = 0;
    end
end

if nargin < 9
    tracefile = [];
end
if nargin < 8
    alg=1;
    params{alg}.levels = 10;
    params{alg}.sol_method = 'MG';
    
    % Top level Acceleration: This can be 'NONE', 'PCG', 'GMRES',
    
    if isSymmetric==0
        params{alg}.TopAcceleration = 'GMRES';
        params{alg}.P_smoothing_omega_val = (5/4)
    else
        params{alg}.TopAcceleration = 'PCG';
        params{alg}.P_smoothing_omega_val = 4/3;
    end
    
    % if NONE is selected in previous option, then "inner" should be 1.
    params{alg}.inner = 20;
    params{alg}.cycle_type = 'V';
    
    % Can be: N1, NN1, ON1, BUp4,
    params{alg}.coarsening_method = 'NN1';
    
    params{alg}.setup_method = 'SA'; % see SetupVcycleT.m
    params{alg}.relax_method = 'GS'; % see GetSmootherInfo.m
    params{alg}.P_smoothing = 'SPAI'; % see smoothProlongationInfo.m
    params{alg}.P_smoothing_omega_type = 'fixed'; % see smoothProlongationInfo.m
    params{alg}.P_epsilonFilter = 0.02; % see Filter.m
    params{alg}.nu1 = 2; % pre relaxations
    params{alg}.nu2 = 2; % post relaxations
    params{alg}.RelParamType = 'fixed'; % see GetSmootherInfo.m
    params{alg}.RelParam = 1;% see GetSmootherInfo.m
    params{alg}.OC = 1.0; % overcorrection parameter
    params{alg}.P_gammaStrong = 0.5; % higher value - larger complexity, better convergence
    params{alg}.oneRelaxOnFineGrid = 1; % do only 1 (pre and post) relaxation on finest grid.
end
if nargin < 7
    TOL = 1e-10;
end
if nargin < 6
    maxit = 200;
end
if nargin < 5
    setup_info = [];
end
setup_info_orig = setup_info;

if nargin < 3 || isempty(x0)
    x0 = zeros(size(b));
    r0 = - b;
end
if isempty(x0)==0
    r0 = AT'*x0-b;
end
rnorm0 = norm(r0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTPUT:

% Generally, x and setup_info are calculated anyway, regardless whether x or
% setup_info are in argout.

% if results is not in argout - no tracing is done...
collecting_results = (nargout == 3); 
% if setting this parameter to be 1 you get all the measurments that were done in the program. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for alg=1:length(params)
    params{alg}.isSymmetric = isSymmetric;
end

if collecting_results
    Tref = timeMATVEC(AT);
    Whole_latex_line = [];
    results = cell(1,length(params));
end

for alg = 1:length(params)
%     disp(num2str(alg));
    
    x = x0;
    r = r0;
    if collecting_results
        results{alg}.rnorms = rnorm0;
        DispToFile(tracefile,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
        DispToFile(tracefile,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
        DispToFile(tracefile,['Initial resdual: ',num2str(rnorm0)]);
        DispToFile(tracefile,['Algorithm: ', num2str(alg),', ',getAlgParamsTitle(params{alg})]);
    end
    if isempty(setup_info_orig)
        [info,trace_info] = doSetup(AT,params{alg});
        setup_info{alg}.info = info;
        setup_info{alg}.trace_info = trace_info;
    else
        info = setup_info{alg}.info;
        trace_info = setup_info{alg}.trace_info;
    end
    if collecting_results
        DispToFile(tracefile,['Cop: ',num2str(trace_info.NNZ/nnz(AT))]);
        DispToFile(tracefile,['Setup took: ',num2str(trace_info.setup_time),' seconds, and ',num2str(round(trace_info.setup_time / Tref)),' MAT-VEC work-units']);
    end
    if collecting_results
        max_nnz = 0;
        mean_nnz = 0;
        Cop_w = 0;
        num_levels = length(info);
        for k = 1:num_levels
            if k==1
                AT_t = AT;
            else
                AT_t = info{k}.AT;
            end
            temp = sum(AT_t~=0);
            max_nnz_t = max(temp);
            max_nnz = max(max_nnz,max_nnz_t);
            mean_nnz_t = mean(temp);
            n = size(AT_t,2);
            if k > 1
                mean_nnz = mean_nnz + mean_nnz_t;
            end
            DispToFile(tracefile,['Level: ',num2str(k),', n = ',num2str(n),', NNZ: ',num2str(nzmax(AT_t)),' Max NNZ per row: ',num2str(max_nnz_t),' Mean: ',num2str(mean_nnz_t)]);
            Cop_w = Cop_w + nzmax(AT_t)*(2^(k-1));
        end
        Cop_w = Cop_w / nzmax(AT);
        mean_nnz = mean_nnz/(num_levels-1);
        
        trace_info.solution_time = 0;
        results{alg}.NNZ = trace_info.NNZ;
        results{alg}.levels = trace_info.levels;
    end
    iter = 0;
    rnorm = rnorm0;
    while iter < maxit
        prev_rnorm = rnorm;
        t = tic;
        [x,r,rnorms,trace_info] = doIteration(AT,x,b,r,params{alg},info,TOL*rnorm0,trace_info);
        t = toc(t);        
        %prev_rnorm = results{alg}.rnorms(end);
        rnorm = rnorms(end); 
        iter = iter + length(rnorms);
        if collecting_results
            trace_info.solution_time = trace_info.solution_time+t;
            results{alg}.rnorms = [results{alg}.rnorms, rnorms];
            gamma = (rnorm/prev_rnorm).^(1/length(rnorms));
            DispToFile(tracefile,['iter: ',num2str(iter),', norm: ',num2str(rnorms(end)),' Conv: ' , num2str(gamma)]);
        end
        
        if rnorm < TOL*rnorm0
            break;
        end
        if (sum(isnan(x))>0)
            error('ERAN: x is nan');
        end
    end
    if collecting_results
        gamma = prod(results{alg}.rnorms(end-8:end)./results{alg}.rnorms(end-9:end-1))^(1/length(results{alg}.rnorms(end-8:end)));
        DispToFile(tracefile,' ');
        DispToFile(tracefile,[params{alg}.setup_method,': Summary for ',getAlgParamsTitle(params{alg})]);
        DispToFile(tracefile,' ');
        Cop = trace_info.NNZ/nnz(AT);
        DispToFile(tracefile,['Iterations: ' ,num2str(iter),' Cop: ',num2str(Cop),', Max NNZ: ',num2str(max_nnz),', Mean: ',num2str(mean_nnz)]);
        gamma_eff = gamma^(1/Cop);
        DispToFile(tracefile,['gamma: ',num2str(gamma), ' gamma_eff ' ,num2str(gamma_eff)]);
        DispToFile(tracefile,['Each MAT-VEC work-unit: ',num2str(Tref),' Seconds.']);
        WUset = round(trace_info.setup_time / Tref);
        WUsol = round(trace_info.solution_time / Tref);
        DispToFile(tracefile,['Setup took: ',num2str(trace_info.setup_time),' seconds, and ',num2str(WUset),' MAT-VEC work-units']);
        DispToFile(tracefile,['Solution took: ',num2str(trace_info.solution_time),' seconds, and ',num2str(WUsol),' MAT-VEC work-units']);
%         latex_line = [num2str(iter),' & ',num2str(Cop,2),' & ', num2str(max_nnz), '& ',num2str(WUset+WUsol)];
%         Whole_latex_line = [Whole_latex_line,' & ',latex_line];
%         DispToFile(tracefile,[latex_line,' \\']);
%         DispToFile(tracefile,' ');
    end
end
% if collecting_results
%     DispToFile(tracefile,' ');
%     DispToFile(tracefile,Whole_latex_line);
%     DispToFile(tracefile,' ');
% end
return;

function [x,r,rnorms,trace_info] = doIteration(AT,x,b,r,params,info,TOL,trace_info)
rnorms = 0;
if strcmp(params.TopAcceleration,'NONE')
    if params.inner > 1
        error('params.inner > 1 && topAcc = NONE');
    end
    x = doInnerIteration(AT,x,b,r,params,info);
    r = AT'*x-b;
    rnorms(end) = norm(r);
%     x = x - doInnerIteration(AT,zeros(length(x),1),r,-r,params,info);%
else
    switch params.TopAcceleration
        case 'GMRES'
            z = zeros(length(x),1);
            prec = @(y)(doInnerIteration(AT,z,y,-y,params,info));
            [x,rnorms,r] = GMRES(AT,x,-r,params.inner,prec,TOL);
        case 'PCG'
            z = zeros(length(x),1);
            prec = @(y)(doInnerIteration(AT,z,y,-y,params,info));
            [x,rnorms,r] = PCG(AT,x,-r,params.inner,prec,TOL);
        case 'RRE'
            X = zeros(length(x),params.inner+1);
            X(:,1) = x;
            for k = 1:params.inner
                X(:,k+1) = doInnerIteration(AT,X(:,k),b,r,params,info);
            end
            x = RRE(X);
            r = AT'*x-b;
    end
end
return;

function x = doInnerIteration(AT,x,b,r,params,info)
switch params.sol_method
    case 'MG'
        x = addVcycleT(x,r,b,params,info,1); 
    case 'Relax'
        x = RelaxT(AT,x,b,r,params.nu1,info,'F');
        x = RelaxT(AT,x,b,r,params.nu1,info,'B');
    otherwise
        error('ERAN: unknown solution method');
end
return;

function [info,trace_info] = doSetup(AT,params)
t = tic;
switch params.sol_method
    case 'Relax'
        trace_info.levels = 1; 
        trace_info.NNZ = nnz(AT);
        info = GetSmootherInfo(params,AT);
    case 'MG'
        [~,info,trace_info] = SetupVcycleT(AT,ones(size(AT,1),1),params,[]);
end
t = toc(t);
trace_info.setup_time = t;
return;
% *******************************************************

function [s] = getAlgParamsTitle(params)
s = [params.TopAcceleration,'(',num2str(params.inner),') acc on '];
if strcmp(params.sol_method,'Relax')
    s = [s, params.sol_method];
else
    s = [s, params.sol_method,',',params.cycle_type,'(',num2str(params.nu1),',',num2str(params.nu2),'); ',params.relax_method ,'(w = ',num2str(params.RelParam),'); ',params.setup_method,'(',params.coarsening_method,')'];
end
return;

function Tref = timeMATVEC(AT)
mm = 100;
y0 = rand(size(AT,2),1);
tic
for k = 1:mm
    r0 = AT*y0;
end
Tref = toc;
Tref = Tref/mm;
return;

