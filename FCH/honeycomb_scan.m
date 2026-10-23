% Honeycomb reference family: hexagonal lattices of light domains with
% spacing d >= 4 (by convention; below that they count as hex spots),
% capped so one hex cell fits in the square sweep box: d*sqrt(3) <= L.
% For each m: polish the lattice at every d, then take the lowest-energy
% one that is stable (even perturbations of the cell).
%
% Settings (override before running): hc_m (0.30:0.02:0.60),
% hc_d (4:0.25:L/sqrt(3)), h (0.098, the resolved sweep spacing), L (4*pi).
% Output: honeycomb_scan.mat with E_honey(m), d_honey(m), and every dE(m, d).

if ~exist('L', 'var'), L = 4*pi; end
if ~exist('hc_m', 'var'), hc_m = 0.30:0.02:0.60; end
if ~exist('hc_d', 'var'), hc_d = 4:0.25:(L/sqrt(3)); end
if ~exist('h', 'var'), h = 0.098; end

Para.epsilon = 0.18;
Para.eta1 = Para.epsilon^2;
Para.eta2 = Para.eta1;
Para.tau = 0;

nm = numel(hc_m);
nd = numel(hc_d);
dE = inf(nm, nd);
stable = false(nm, nd);
E_honey = inf(1, nm);
d_honey = nan(1, nm);
q = parallel.pool.DataQueue;
t_start = tic;
afterEach(q, @(msg) fprintf('[%5.1f min] %s\n', toc(t_start)/60, msg));
parfor a = 1:nm
    m = hc_m(a);
    e = inf(1, nd);
    states = cell(1, nd);
    for j = 1:nd
        [u, G, info, kind] = SSAV_FCH_helpers.hex_lattice_state(m, hc_d(j), 'light', Para, h);
        if strcmp(kind, 'lattice')
            e(j) = info.dE;
            states{j} = struct('u', u, 'G', G);
        end
    end
    % lowest stable one: check stability in order of energy
    st = false(1, nd);
    [~, ord] = sort(e);
    for j = ord(isfinite(e(ord)))
        [~, ~, si] = stability_FCH(states{j}.u, states{j}.G, Para, struct('even_only', true, 'verbose', false));
        st(j) = si.stable;
        if st(j), break, end
    end
    dE(a, :) = e;
    stable(a, :) = st;
    es = e;  es(~st) = inf;
    [E_honey(a), jb] = min(es);
    if isfinite(E_honey(a)), d_honey(a) = hc_d(jb); end
    send(q, sprintf('m = %.2f: best stable honeycomb dE %+.6f at d = %.2f', m, E_honey(a), d_honey(a)));
end
save('honeycomb_scan.mat', 'hc_m', 'hc_d', 'dE', 'stable', 'E_honey', 'd_honey', 'h', 'L', 'Para');
fprintf('Saved honeycomb_scan.mat\n');
