function Seqs = Simulation_Thinning_HP_Safe(para, options)
%SIMULATION_THINNING_HP_SAFE Simulate options.N non-empty event sequences.
% Wrapper around the toolkit's Simulation_Thinning_HP (Ogata thinning):
% sequences without events are discarded and regenerated until options.N
% non-empty sequences are available (or 2*options.N attempts are made).

max_attempts = options.N * 2;
Seqs = [];
attempts = 0;

while length(Seqs) < options.N && attempts < max_attempts
    opts_batch = options;
    opts_batch.N = options.N - length(Seqs);

    Seqs_batch = Simulation_Thinning_HP(para, opts_batch);

    % keep sequences with at least one event
    valid_seqs = arrayfun(@(s) ~isempty(s.Time), Seqs_batch);
    Seqs = [Seqs, Seqs_batch(valid_seqs)]; %#ok<AGROW>

    empty_count = sum(~valid_seqs);
    if empty_count > 0
        fprintf('  Removed %d empty sequences, regenerating...\n', empty_count);
    end

    attempts = attempts + opts_batch.N;
end

if length(Seqs) < options.N
    warning('Could not generate %d valid sequences after %d attempts', options.N, attempts);
end

fprintf('Generated %d valid sequences\n', length(Seqs));
end
