function [Qcg,omega_cg] = smoothProlongationInfo(params,AT,P0)
cgParams.relax_method = params.P_smoothing;
cgParams.RelParamType =  params.P_smoothing_omega_type;
cgParams.RelParam = params.P_smoothing_omega_val;
cgInfo = GetSmootherInfo(cgParams,AT);
Qcg = cgInfo.Q;

rho = max(Qcg./GetDiagPreconditioner(AT,'1'));

omega_cg = cgInfo.omega/rho;

return;