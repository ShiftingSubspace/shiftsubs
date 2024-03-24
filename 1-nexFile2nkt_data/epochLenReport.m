function epochLength = epochLenReport(MT_cond, epochStatDesp, mkIdDict, reportMode)
% report the length of each behavior epoch
if nargin < 4 
    reportMode = 'v';
end

MT_cond = [cell2mat(MT_cond{1}); cell2mat(MT_cond{2}); cell2mat(MT_cond{3}); cell2mat(MT_cond{4})];
epochs = fieldnames(epochStatDesp);
epochLength = nan(1, length(epochs));
for ep =1:length(epochs)
    currLeft = epochStatDesp.(epochs{ep}).('left');
    currRight = epochStatDesp.(epochs{ep}).('right');
    currLeftIndex = mkIdDict.(currLeft);
    currRightIndex = mkIdDict.(currRight);
    t_diff = MT_cond(:, currRightIndex) - MT_cond(:, currLeftIndex);
    epochLength(ep) = median(t_diff);

    if isequal(reportMode, 'v')
        fprintf("%s - %s: %.2f\n", currLeft, currRight, epochLength(ep));
    end
end

