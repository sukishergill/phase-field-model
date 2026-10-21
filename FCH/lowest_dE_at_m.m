function [e, which] = lowest_dE_at_m(branches, mq, stable_only)
% Lowest dE over a set of continuation branches at each mass in mq.
%
%   branches     cell array of continuation_FCH branch structs (any shape)
%   mq           vector of masses
%   stable_only  only use stretches with n_neg = 0 at both ends (true)
%
%   e      lowest dE at each mq (inf where no branch covers that m),
%          interpolated linearly in m between neighbouring branch points
%   which  linear index into branches of the branch attaining it (NaN if none)

if nargin < 3, stable_only = true; end
e = inf(size(mq));
which = nan(size(mq));
for a = 1:numel(branches)
    b = branches{a};
    m1 = b.m(1:end-1);  m2 = b.m(2:end);
    e1 = b.dE(1:end-1); e2 = b.dE(2:end);
    ok = m1 ~= m2;
    if stable_only
        ok = ok & b.n_neg(1:end-1) == 0 & b.n_neg(2:end) == 0;
    end
    m1 = m1(ok); m2 = m2(ok); e1 = e1(ok); e2 = e2(ok);
    for q = 1:numel(mq)
        in = (m1 - mq(q)) .* (m2 - mq(q)) <= 0;
        if ~any(in), continue, end
        t = (mq(q) - m1(in)) ./ (m2(in) - m1(in));
        ek = min(e1(in) + t .* (e2(in) - e1(in)));
        if ek < e(q)
            e(q) = ek;
            which(q) = a;
        end
    end
end
end
