function [branches, files] = load_sweep_branches(pattern)
% Stack the 'other' (m-sweep) branches of all continue_m_sweep outputs
% (each an nm x 2 cell: starting state x direction) into one cell array,
% so several sweeps are compared together.
%
%   [branches, files] = load_sweep_branches()          % all *_tofold_branches.mat
%   [branches, files] = load_sweep_branches(pattern)   % e.g. 'ex5_m_sweep_*_branches.mat'

if nargin < 1, pattern = 'ex5_m_sweep_*_tofold_branches.mat'; end
f = dir(pattern);
files = {f.name};
branches = cell(0, 2);
for k = 1:numel(files)
    S = load(files{k}, 'branches');
    branches = [branches; S.branches]; %#ok<AGROW>
end
% order by starting mass
m0 = cellfun(@(b) b.m(1), branches(:, 1));
[~, ord] = sort(m0);
branches = branches(ord, :);
end
