function ExampleIdx = SelectExample(Params, Examples)
% function ExampleIdx = SelectExample(Params, Examples)

% (c) XXXX

while true
  if ~isfield(Params,'ExampleName')
    fprintf('\n Examples:\n');
    for k = 1:length(Examples)
      fprintf('\n [%2d] %s',k,Examples{k}.sysInfo.name);
    end
    fprintf('\n\n');    
    ExampleIdx = input('Pick an example to run:           ');
    try
      fprintf('\nRunning %s\n',Examples{ExampleIdx}.sysInfo.name);
      break;
    catch
    end
  else
    ExampleIdx  = find(cellfun(@(x) strcmp(x.sysInfo.name, Params.ExampleName), Examples));
    if isempty(ExampleIdx)
      error('SOD_Utils:SelectExample:exception', 'ExampleName in Params has no matching entry in all the pre-set examples!!');
    else
      fprintf('\nRunning %s\n', Examples{ExampleIdx}.sysInfo.name);
      break;
    end
  end
end

return