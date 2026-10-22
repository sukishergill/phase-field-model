function Rf = load_reference_patterns(files)
% Load and merge reference pattern branches (reference_patterns_sweep.m).
% Every file holds valid steady states, so runs with different settings
% (e.g. reference_patterns_v1.mat, before the uniform-crossing stop, and
% the current reference_patterns.mat) are combined: the rows of stripes /
% spots are stacked, with n_stripes / spot_d labels to match.
%
%   Rf = load_reference_patterns()          % all reference_patterns*.mat
%   Rf = load_reference_patterns({'a.mat'}) % given files
%
% Hex-cell branches can jump onto other patterns that fit the cell -- with
% ds_max = 0.2 one d = 3.14 branch ended up on two straight stripes, which
% then counted as "hex spots". Points of a spot branch that are stripes
% (u constant along one axis to within 5% of its variation along the other)
% get dE = NaN, so they drop out of comparisons and plots; Rf.n_dropped
% counts them.

if nargin < 1
    f = dir('reference_patterns*.mat');
    files = {f.name};
end
Rf = struct('stripes', {cell(0, 2)}, 'n_stripes', [], 'spots', {cell(0, 2)}, ...
    'spot_d', [], 'files', {files}, 'n_dropped', 0);
for k = 1:numel(files)
    S = load(files{k}, 'stripes', 'n_stripes', 'spots', 'spot_d');
    Rf.stripes = [Rf.stripes; S.stripes];
    Rf.n_stripes = [Rf.n_stripes, S.n_stripes];
    for j = 1:numel(S.spots)
        b = S.spots{j};
        for k = 1:numel(b.m)
            u = b.u{k};
            vx = mean(std(u, 0, 2));    % variation along x
            vy = mean(std(u, 0, 1));    % variation along y
            if min(vx, vy) < 0.05*max(vx, vy)
                b.dE(k) = NaN;          % a stripe, not a spot lattice
                Rf.n_dropped = Rf.n_dropped + 1;
            end
        end
        S.spots{j} = b;
    end
    Rf.spots = [Rf.spots; S.spots];
    Rf.spot_d = [Rf.spot_d, S.spot_d];
end
end
