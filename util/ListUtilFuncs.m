function ListUtilFuncs (prefix)

if nargin==0, prefix = '..'; end

anlypath = [prefix '/utils/anly/'];
anlyfiles = dir([anlypath '*.m']);

plotpath = [prefix '/utils/plot/'];
plotfiles = dir([plotpath '*.m']);


fprintf('\n')
fprintf('Functions in %s:\n', anlypath);
fprintf('%s\n', anlyfiles.name);


fprintf('\n')
fprintf('Functions in %s:\n', plotpath);
fprintf('%s\n', plotfiles.name);
fprintf('\n');

end