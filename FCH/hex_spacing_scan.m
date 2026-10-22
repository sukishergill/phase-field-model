% Hexagonal spot lattice energy against lattice spacing d at fixed masses,
% to check whether any spacing beats a given "other" state.
%
% Each lattice is computed on one exact hex cell (width d, height d*sqrt(3),
% two spots, kept identical by the half-cell shift symmetry), started from
% spots at the target m and polished there (newton_polish_FCH, even). Both
% polarities: dark spots (u ~ -1) in u ~ +1, and light spots in u ~ -1.
% A result is flagged if it is not a spot lattice: stripes (constant along
% an axis) or uniform. Stability is within even perturbations of the cell.
%
% Settings (override before running):
%   scan_m      masses                        ([0.026 0.118])
%   scan_d      spacings                      (1.2:0.05:3.2)
%   other_dE    "other" dE at each mass, for the plot (resolution_check
%               values for the m-sweep states at 0.026 and 0.118)
%   h           target grid spacing           (0.08)

if ~exist('scan_m', 'var'), scan_m = [0.026 0.118]; end
if ~exist('scan_d', 'var'), scan_d = 1.2:0.05:3.2; end
if ~exist('other_dE', 'var'), other_dE = [-0.0002125736, -0.0058335696]; end
if ~exist('h', 'var'), h = 0.08; end

Para.epsilon = 0.18;
Para.eta1 = Para.epsilon^2;
Para.eta2 = Para.eta1;
Para.tau = 0;
ep = Para.epsilon;
layer = @(s) (1 - tanh(s / (sqrt(2)*ep))) / 2;   % ~1 inside (s < 0)
half_shift = @(v) (v + circshift(v, size(v)/2)) / 2;
popts = struct('even', true, 'symmetry', half_shift, 'verbose', false, 'maxit', 3000);

nm = numel(scan_m);  nd = numel(scan_d);
res = struct('dE', nan(nm, nd, 2), 'n_neg', nan(nm, nd, 2), 'ok', false(nm, nd, 2), ...
    'kind', {strings(nm, nd, 2)}, 'converged', false(nm, nd, 2));
pol_name = {'dark spots', 'light spots'};

for a = 1:nm
    m = scan_m(a);
    for j = 1:nd
        d = scan_d(j);
        Ly = d*sqrt(3);
        Nx = 2*ceil(d / (2*h));
        Ny = 2*ceil(Ly / (2*h));
        G = SSAV_FCH_helpers.generate_Grid([d, Ly], [Nx, Ny], 2);
        r1 = hypot(G.xx, G.yy);
        ddx = G.xx - d/2;   ddx = ddx - d*round(ddx/d);
        ddy = G.yy - Ly/2;  ddy = ddy - Ly*round(ddy/Ly);
        r2 = hypot(ddx, ddy);
        for p = 1:2
            % minority-phase area fraction at this m sets the initial radius
            f = (1 - (3 - 2*p)*m) / 2;                      % dark: (1-m)/2, light: (1+m)/2
            R = sqrt(f * d * Ly / (2*pi));
            spots = layer(r1 - R) + layer(r2 - R);
            if p == 1, u = 1 - 2*spots; else, u = -1 + 2*spots; end
            u = u - mean(u(:)) + m;
            [u, info] = newton_polish_FCH(u, G, Para, popts);
            vx = mean(std(u, 0, 2));  vy = mean(std(u, 0, 1));
            if info.collapsed
                kind = "uniform";
            elseif min(vx, vy) < 0.05*max(vx, vy)
                kind = "stripes";
            else
                kind = "spots";
            end
            [~, ~, si] = stability_FCH(u, G, Para, struct('even_only', true, 'verbose', false));
            res.dE(a, j, p) = info.dE;
            res.n_neg(a, j, p) = si.n_neg;
            res.kind(a, j, p) = kind;
            res.converged(a, j, p) = info.converged;
            res.ok(a, j, p) = info.converged && kind == "spots";
        end
    end
    % report the best spot lattice at this m
    fprintf('m = %.3f:\n', m);
    for p = 1:2
        e = res.dE(a, :, p);  e(~res.ok(a, :, p)) = NaN;
        [eb, jb] = min(e);
        es = e;  es(res.n_neg(a, :, p) > 0) = NaN;
        [esb, jsb] = min(es);
        fprintf('  %-11s: best dE %+.7f at d = %.2f (n_neg %d); best stable %+.7f at d = %.2f | %d spots, %d stripes, %d uniform, %d not converged\n', ...
            pol_name{p}, eb, scan_d(jb), res.n_neg(a, jb, p), esb, scan_d(jsb), ...
            sum(res.kind(a, :, p) == "spots"), sum(res.kind(a, :, p) == "stripes"), ...
            sum(res.kind(a, :, p) == "uniform"), sum(~res.converged(a, :, p)));
    end
    e_all = res.dE(a, :, :);  e_all(~res.ok(a, :, :) | res.n_neg(a, :, :) > 0) = NaN;
    fprintf('  "other" dE %+.7f  ->  other - best stable hex = %+.3e  (%s)\n', other_dE(a), ...
        other_dE(a) - min(e_all(:)), ternary(other_dE(a) < min(e_all(:)), 'other lower', 'a hex spacing is lower'));
end
save('hex_spacing_scan.mat', 'res', 'scan_m', 'scan_d', 'other_dE', 'h', 'Para');

%% plot dE against d
f = figure('Position', [100 100 1100 420]);
tl = tiledlayout(f, 1, nm, 'TileSpacing', 'compact');
cols = [0.75 0 0.75; 0.2 0.5 0.9];
for a = 1:nm
    ax = nexttile(tl); hold(ax, 'on');
    hp = gobjects(2, 1);
    for p = 1:2
        e = res.dE(a, :, p);  e(~res.ok(a, :, p)) = NaN;
        st = e;  st(res.n_neg(a, :, p) > 0) = NaN;     % stable: filled
        hp(p) = plot(ax, scan_d, e, '-', 'Color', cols(p, :), 'LineWidth', 1);
        plot(ax, scan_d, st, 'o', 'Color', cols(p, :), 'MarkerFaceColor', cols(p, :), 'MarkerSize', 4);
    end
    ho = yline(ax, other_dE(a), 'k--', 'LineWidth', 1.5);
    legend(ax, [hp; ho], [strcat("hex, ", pol_name), "other (m-sweep state)"], 'Location', 'best');
    xlabel(ax, 'lattice spacing d');
    ylabel(ax, 'E - E_{uniform}');
    title(ax, sprintf('m = %.3f (filled: stable)', scan_m(a)));
    grid(ax, 'on'); box(ax, 'on');
end
exportgraphics(f, 'hex_spacing_scan.png', 'Resolution', 110);
fprintf('Saved hex_spacing_scan.mat and hex_spacing_scan.png\n');

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
