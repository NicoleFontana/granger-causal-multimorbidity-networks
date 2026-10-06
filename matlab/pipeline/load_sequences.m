function [Seqs, mapping] = load_sequences(analysis_dir)
%LOAD_SEQUENCES Read sequences.csv (patient_id, Time, Mark) and
% event_mapping.csv written by R/01_build_sequences.R, and convert them to
% the Hawkes-Process-Toolkit sequence structure (Time, Mark, Start, Stop,
% Feature). Time is in months since age 40; the observation window of each
% individual ends at the last recorded event.

T = readtable(fullfile(analysis_dir, 'sequences.csv'));
mapping = readtable(fullfile(analysis_dir, 'event_mapping.csv'), 'TextType', 'string');

T = sortrows(T, {'patient_id', 'Time', 'Mark'});
[ids, ~, g] = unique(T.patient_id);
time_by = splitapply(@(x) {x(:)'}, T.Time, g);
mark_by = splitapply(@(x) {x(:)'}, T.Mark, g);

n = numel(ids);
Seqs(n) = struct('Time', [], 'Mark', [], 'Start', 0, 'Stop', [], 'Feature', []);
for i = 1:n
    Seqs(i).Time    = time_by{i};
    Seqs(i).Mark    = mark_by{i};
    Seqs(i).Start   = 0;
    Seqs(i).Stop    = max(time_by{i});
    Seqs(i).Feature = [];
end
end
