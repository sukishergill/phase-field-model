function [u, hist, pinfo] = descend_FCH(u, Grid, Para, opts)
% From a polished FCH steady state, descend to a stable one: while the
% state has a negative eigenvalue (stability_FCH), perturb it along the most
% negative mode in both directions, re-polish each with PTC
% (newton_polish_FCH), and keep the lower-energy converged result.
%
% PTC only accepts energy-decreasing steps, so it cannot slide back to the
% saddle once pushed off it: E(saddle + delta*d) < E(saddle) when d has
% negative curvature.
%
% Inputs
%   u, Grid, Para  as for newton_polish_FCH (u should be polished)
%   opts (optional)
%     delta        RMS size of the push along the unstable mode (0.02)
%     max_rounds   max saddle-escape rounds                     (10)
%     polish_opts  options passed to newton_polish_FCH  (maxit 3000, quiet)
%     stab_opts    options passed to stability_FCH      (quiet)
%     verbose      print progress                               (true)
%
% Outputs
%   u     final state (stable unless max_rounds ran out or a polish failed)
%   hist  struct array per round: E, n_neg, lam (smallest eigenvalues)
%   pinfo newton_polish_FCH info for the final state

if nargin < 4, opts = struct(); end
defaults = struct('delta', 0.02, 'max_rounds', 10, ...
    'polish_opts', struct('verbose', false, 'maxit', 3000), ...
    'stab_opts', struct('verbose', false), 'verbose', true);
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}), opts.(names{i}) = defaults.(names{i}); end
end

[u, pinfo] = newton_polish_FCH(u, Grid, Para, opts.polish_opts);
hist = struct('E', {}, 'n_neg', {}, 'lam', {});

for round = 1:opts.max_rounds + 1
    [lam, modes, sinfo] = stability_FCH(u, Grid, Para, opts.stab_opts);
    hist(end+1) = struct('E', pinfo.E, 'n_neg', sinfo.n_neg, 'lam', lam); %#ok<AGROW>
    if opts.verbose
        fprintf('  round %d: E = %+.10f, %d negative eigenvalue(s), smallest %s\n', ...
            round - 1, pinfo.E, sinfo.n_neg, mat2str(lam(1:min(3, end))', 3));
    end
    if sinfo.stable || round > opts.max_rounds
        break
    end

    % push along the most negative mode, both signs; keep the lower energy
    d = modes(:, :, 1);
    best = [];
    for s = [1 -1]
        [u_try, info_try] = newton_polish_FCH(u + s*opts.delta*d, Grid, Para, opts.polish_opts);
        if info_try.converged && info_try.E < pinfo.E && (isempty(best) || info_try.E < best.info.E)
            best = struct('u', u_try, 'info', info_try);
        end
    end
    if isempty(best)
        if opts.verbose
            fprintf('  no converged lower-energy state found from either push -- stopping\n');
        end
        break
    end
    u = best.u;
    pinfo = best.info;
end
end
