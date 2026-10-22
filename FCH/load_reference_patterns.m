function Rf = load_reference_patterns(files)
% Load and merge reference pattern branches (reference_patterns_sweep.m).
% Every file holds valid steady states, so runs with different settings
% (e.g. reference_patterns_v1.mat, before the uniform-crossing stop, and
% the current reference_patterns.mat) are combined: the rows of stripes /
% spots are stacked, with n_stripes / spot_d labels to match.
%
%   Rf = load_reference_patterns()          % all reference_patterns*.mat
%   Rf = load_reference_patterns({'a.mat'}) % given files

if nargin < 1
    f = dir('reference_patterns*.mat');
    files = {f.name};
end
Rf = struct('stripes', {cell(0, 2)}, 'n_stripes', [], 'spots', {cell(0, 2)}, ...
    'spot_d', [], 'files', {files});
for k = 1:numel(files)
    S = load(files{k}, 'stripes', 'n_stripes', 'spots', 'spot_d');
    Rf.stripes = [Rf.stripes; S.stripes];
    Rf.n_stripes = [Rf.n_stripes, S.n_stripes];
    Rf.spots = [Rf.spots; S.spots];
    Rf.spot_d = [Rf.spot_d, S.spot_d];
end
end
