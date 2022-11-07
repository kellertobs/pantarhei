function [params, flds] = ListRunProperties (folder, RunIDVec)

Nf = length(RunIDVec);

if ~iscell(folder), folder = {folder}; end
if length(folder)<Nf, folder = repmat(folder, Nf, 1); end

for fi = 1:Nf
    fp = GetOutputMatFiles(folder{fi}, RunIDVec{fi});
    params(:,fi) = struct2cell(load(fp));
end

flds = fieldnames(load(fp));

array2table(params, 'VariableNames', RunIDVec', 'RowNames', flds)

end