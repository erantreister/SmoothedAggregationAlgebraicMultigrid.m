function PT_filtered = ProlongationFilter(PT,xc,epsilon)
% no X:
% epsilon = 0;
% PT = [1   ,0   ,0,0,1    ,1   ,3   ,4;
%      0.1 ,1e-4 ,3,1,1e-2 ,0.5 ,1e-3,0;
%      0   ,3    ,0,1,1e-4 ,1e-5,0   ,1e-5];
ABSPT = abs(PT);
maxPT = max(ABSPT)';
PT_normalized = ABSPT*toSparseDiag(1./maxPT);
I_chosen = PT_normalized > epsilon; clear ABSPT; clear PT_normalized;

PT_chosen = PT.*I_chosen;
sumPT_chosen = sum(PT_chosen)';
weights_chosen = PT_chosen*toSparseDiag(1./sumPT_chosen);
PT_collapsed = weights_chosen*toSparseDiag((sum(PT) - sum(PT_chosen))');
PT_filtered = PT_chosen + PT_collapsed;
% norm(sum(PT_filtered)-sum(PT))

end

