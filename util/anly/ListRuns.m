function fname = ListRuns (folder)

if nargin==0
    sweeppath = '../pantarhei/out/sweep/';
else
    sweeppath = folder;
end

% get directories in folder path
sweepfiles = dir(sweeppath);

% remove system directories
sweepfiles(contains({sweepfiles.name}, '.')) = [];

if nargout==0
    % just print out the folders
    fprintf('\n')
    fprintf('Folders in %s:\n', sweeppath);
    fprintf('%s\n', sweepfiles.name);
    fprintf('\n')
else
    % output to variable
    fname = {sweepfiles.name}';
end

end