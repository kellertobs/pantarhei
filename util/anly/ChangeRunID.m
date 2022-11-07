function [newpath, newfnames] = ChangeRunID (folder, oldID, newID, deloldID)
% changeRunID (folder, oldID, newID, deloldID)
% 
% example: changeRunID('../out/', 'test', 'olv10_plg10_bas80', 0);
% 
% helper function to change the RunID of a set of simulations.
% 
% INPUTS
% folder    where the run is stored, usually '../out/'
% oldID     old simulation name
% newID     new simulation name
% deloldID  whether you want to delete the old directory. default = no
%           Be Careful!!
% 
% OUTPUTS (optional)
% newpath   path of new directory of simulation
% newfnames new file names in the new directory
% 
% YQW, 6 May 2021


% define file paths
oldpath = [folder, oldID '/'];
newpath = [folder, newID '/'];

newfnames = {};

% make new folder
if ~exist(newpath,'dir')
    mkdir(newpath);
else
    error('!! %s already exists. Choose another RunID !!\n\n', newpath);
end
    

% get files in old directory
tmp       = dir(oldpath);
oldfnames = {tmp.name}';
oldfnames(strcmp(oldfnames, '.')) = [];
oldfnames(strcmp(oldfnames, '..')) = [];
%oldfnames = oldfnames(contains(oldfnames, '.pdf') | contains(oldfnames, '.mat') | contains(oldfnames, '.mp4'));

% replace oldID with newID in filenames
newfnames = strrep(oldfnames, oldID, newID);

% now copy files
for fi = 1:length(oldfnames)
    copyfile(strcat(oldpath, oldfnames{fi}), strcat(newpath, newfnames{fi}));
end


if nargin<4, deloldID = 0; end

% flag for deleting the old directory
if (deloldID)
    rmdir([folder, oldID], 's');
end


end