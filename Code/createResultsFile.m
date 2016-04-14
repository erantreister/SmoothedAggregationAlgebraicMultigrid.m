function file = createResultsFile(filename)
filename = ['./results/',filename,'-',num2str(round(rand()*10000))];
a = [];
save([filename,'.txt'],'a','-ASCII');
file = fopen([filename,'.txt'],'at');
return;


% function file = createResultsFile(matName,unknowns_size)
% filename = ['./results/',matName,'-',num2str(unknowns_size),'-',num2str(round(rand()*1000))];
% a = [];
% save([filename,'.txt'],'a','-ASCII');
% file = fopen([filename,'.txt'],'at');
% return;