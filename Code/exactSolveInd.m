function [res] = exactSolveInd(n,lev,params)
res = (n <= 200 || params.levels <= lev);
return;