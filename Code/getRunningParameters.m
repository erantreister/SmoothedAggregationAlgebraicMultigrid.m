function [params] = getRunningParameters(MatName)
params = cell(1);
nu1 = 1;
nu2 = 1;

epsilonFilter = 0.02;
gammaStrong = 0.5;
Psmoothing_omega_val = 4/3;
RelParam = 1;



alg = 0;
% alg=alg+1;
% params{alg}.levels = 1;
% params{alg}.TopAcceleration = 'GMRES';
% params{alg}.inner = 100;
% params{alg}.sol_method = 'Relax';
% params{alg}.relax_method = 'GS';
% params{alg}.nu = 10; 
% params{alg}.RelParamType = 'fixed';
% params{alg}.RelParam = 1;


alg=alg+1;
params{alg}.levels = 15;
params{alg}.sol_method = 'MG';      % {Relax}
params{alg}.TopAcceleration = 'PCG';% {GMRES, PCG, NULL}
params{alg}.inner = 50;             % {if NULL: should be 1}
params{alg}.cycle_type = 'V';       % {V,W,F} 
params{alg}.coarsening_method = 'NN1';% see getAggregation.m
params{alg}.setup_method = 'AGG';
params{alg}.relax_method = 'GS';
params{alg}.nu1 = nu1;
params{alg}.nu2 = nu2;
params{alg}.RelParamType = 'fixed';
params{alg}.RelParam = RelParam;
params{alg}.OC = 'rec_minres';
params{alg}.P_gammaStrong = gammaStrong;
params{alg}.oneRelaxOnFineGrid = 0;

% alg=alg+1;
% params{alg} = params{alg-1};
params{alg}.setup_method = 'SpSA';
params{alg}.P_smoothing = 'SPAI';
params{alg}.P_smoothing_omega_type = 'fixed';
params{alg}.P_smoothing_omega_val = Psmoothing_omega_val;
params{alg}.P_epsilonFilter = epsilonFilter;
params{alg}.OC = 1.0;

% alg=alg+1;
% params{alg} = params{alg-1};
params{alg}.setup_method = 'SA';
params{alg}.P_smoothing = 'SPAI';
params{alg}.P_smoothing_omega_type = 'fixed';
params{alg}.P_smoothing_omega_val = Psmoothing_omega_val;
params{alg}.P_epsilonFilter = epsilonFilter;

return;