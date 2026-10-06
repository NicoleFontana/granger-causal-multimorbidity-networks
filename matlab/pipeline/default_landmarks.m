function landmark = default_landmarks(L, Tmax)
%DEFAULT_LANDMARKS Centres (months) of the L Gaussian basis functions.
% Denser at short lags: 40% of the landmarks equally spaced in 0-12 months,
% 30% log-spaced in 12-60 months and the remainder log-spaced in 60-120
% months (capped at Tmax). For L = 10 this gives 4 + 3 + 3 landmarks.

n_short  = ceil(L * 0.4);
n_medium = ceil(L * 0.3);
n_long   = L - n_short - n_medium;

short_lm  = linspace(0, 12, n_short);
medium_lm = logspace(log10(12.1), log10(60), n_medium);
long_lm   = logspace(log10(60.1), log10(min(120, Tmax)), n_long);

landmark = sort(unique([short_lm, medium_lm, long_lm]));
end
