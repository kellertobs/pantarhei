function [p] = ExtractParams (folder, RunIDVec, params)
% 
% v = ExtractParams (folder, RunIDVec, vars)
% 
% example: v = ExtractParams('../out/', {'test1','test2'}, {'A','B','C'})
% 
% Extract simulation parameters from the par.mat file
% 
% INPUTS
% folder  	folder name where the simulations are stored
% RunIDVec 	vector of simulation names [Nruns x 1]
% params    cell vector of parameter names [Npars x 1]
% 
% OUTPUTS
% p         cell vector of run parameters {Nruns x Npars}
% 
% YQW, 19 Jan 2021

% get vector lengths
Nruns = length(RunIDVec);
Npars = length(params);

% initialise
p = cell(Nruns, Npars);

% assign parameters
for ri = 1:Nruns
    fp      = GetOutputMatFiles(folder, RunIDVec{ri});
    p(ri,:) = struct2cell(load(fp, params{:}));
end

end