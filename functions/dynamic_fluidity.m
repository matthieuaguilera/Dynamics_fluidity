% function: Fluidity = dynamic_fluidity(TS)
%
% Computes the dynamic fluidity of a multivariate time series as introduced in 
% Aguilera et al. (2025). This metric quantifies local irregularity in state space 
% transitions based on extreme value theory (EVT). Time points with locally extreme 
% divergence patterns—i.e., whose distance to surrounding states is highly uneven—
% yield higher fluidity values.
%
% INPUT:
%   TS        : Time series matrix. Can be in Time × Variables or Variables × Time format.
%               The function automatically detects the orientation and assumes the 
%               smallest dimension corresponds to the number of variables.
%
% OUTPUT:
%   Fluidity  : Row vector of dynamic fluidity values, one per time point.
%
% REFERENCES:
%   Faranda, D., Messori, G., & Yiou, P. (2017). Dynamical proxies of North Atlantic 
%   predictability and extremes. Scientific Reports, 7, 41278. 
%   https://doi.org/10.1038/srep41278
%
%   Aguilera, M., Mathis, C., Herbeaux, K., Isik, A., Faranda, D., Battaglia, D., 
%   Goutagny, R. (2025). 40 Hz light stimulation restores early brain dynamics 
%   alterations and associative memory in Alzheimer’s disease model mice. 
%   Imaging Neuroscience, 3, IMAG.a.70. https://doi.org/10.1162/IMAG.a.70
%
% Author: Matthieu Aguilera — FunSy team, LNCA (Strasbourg), 2025

function Fluidity = dynamic_fluidity(TS)

% --- Parameters ---
quanti = 0.98;   % Quantile level for EVT-based estimation
STEP = 1;        % Subsampling step size (keep = 1 for full resolution)

% --- Handle input orientation ---
[rows, cols] = size(TS);
if cols > rows  % Assume input is Variables × Time
    TS = TS';   % Transpose to Time × Variables
end

N = size(TS, 1);     % Number of time points
Fluidity = nan(1, N);  % Preallocate

% --- Main loop ---
for j = 1:STEP:N
    others = setdiff(1:STEP:N, j);
    distance = pdist2(TS(j, :), TS(others, :));
    logdista = -log(distance);
    Fluidity(j) = extremal_Sueveges(logdista, quanti);
end

end

% --- EVT estimator for local irregularity ---
% Estimate the extremal index using Sueveges' maximum likelihood method
function theta = extremal_Sueveges(Y, p)

u = quantile(Y, p);
q = 1 - p;
Li = find(Y > u);
Ti = diff(Li);
Si = Ti - 1;

Nc = sum(Si > 0);     % number of non-consecutive exceedances
N = length(Ti);       % total exceedances

S = sum(q .* Si);
theta = (S + N + Nc - sqrt((S + N + Nc)^2 - 8 * Nc * S)) / (2 * S);
end