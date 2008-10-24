function change_params_file(oldparamfn, newparamfn, pnames, pvalues)
%function change_params_file(oldparamfn, newparamfn, pnames, pvalues)
%
% function for modifying an ascii-parameter file by setting various
% parameters. The oldparamfn file is read and the newparamfn is generated
% or overwritten without check. All lines setting any of the 
% parameters listed in the cell array pnames (with string entries) 
% are removed from the input-param-file and are replaced with a line
% setting the new values given in the cell array pvalues (containing strings)
%
% a sample call would be:
%
% change_params_file('potentialparams.dat','newpotentialparams.dat',...
%                     {'ElectronicDirichletValue','Grid'},...
%                     {'0.0', 'mypathtomacrofiles/othermacro.hexa'});
  
% Bernard Haasdonk 23.2.2007

% open input file
oldp = textread('potentialparam.dat','%[^\n]','commentstyle','matlab')
    
% open output file
  
fid = fopen(newparamfn,'wt');

% for all lines

for nline = 1:length(oldp)
  %   check, whether parameter is to be changed
  p = oldp{nline};
  i = findstr(p,':');
  pname = p(1:(min(i)-1));
  %disp(['read paramname: ',pname]);  
  %   either write line unmodified or filter line (do nothing)
  if ~ismember({pname},pnames)
    fprintf(fid,'%s\n',p);
  else
    fprintf(fid,'%%%%%% The following parameter is replaced by change_params_file.m\n');
    fprintf(fid,'%% %s\n',p);
  end;
  % keyboard;
end;

% for all parameters to be set, insert line into file  
fprintf(fid,'%%%%%% The following entries are generated by change_params_file.m\n');
for npar = 1:length(pnames)
  fprintf(fid,'%s: %s \n',pnames{npar}, pvalues{npar});
end;

% close files
fclose(fid);
  
  
 