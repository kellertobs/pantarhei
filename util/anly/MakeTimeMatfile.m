function [t, tfname] = MakeTimeMatfile (folder, RunID)
% 
% [t, matfilename] = MakeTimeMatfile (folder, RunID)
%
% example:
% [t, matfilename] = MakeTimeMatfile ('../out/', 'olv10_plg10_bas80');
% 
% makes a mat file with the list of time points of the simulation run. 
% To be used at the end of the simulation, to collect the time points so
% that you can query time-specific mat files rather than load all the
% output mat files (takes a long time)
% 
% 
% INPUTS
% folder  	folder name where the simulation is stored
% RunID 	name of simulation
%
% OUTPUTS
% t         time of simulation [Nf x 1]
% tfilename name of the output mat file
% 
% YQW, 22 Aug 2022
% 

% get files from simulations
[~, fn] = GetOutputMatFiles(folder, RunID);
Nf      = length(fn);          % number of files

% initialise outputs
tfname = [folder, RunID, '/', RunID '_tvec.mat'];
t      = zeros(Nf,1);

for fi = 1:Nf
    tmp   = load(fn{fi}, 'time');
    t(fi) = tmp.time;
end

save(tfname, 't');

end