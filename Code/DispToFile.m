function DispToFile(file,s)    
disp(s);
if isempty(file)==0,
    fprintf(file, '%s\n',s);
end
return;
