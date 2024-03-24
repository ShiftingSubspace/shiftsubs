function [binnedSpikesOut, conf] = dpcabin_formatSpikes(spk, conf)
%% spikes into bin with input config specified in conf. If not specified, will use the default conf. 
if ~exist('conf','var')
    conf.prior = 1.0;
    conf.after = 1.0;
    conf.tstep = 1/1000;
end
if iscell(spk)

        binnedSpikesOut = cellfun(@(x) dpcabin_formatSpikes(x, conf),spk,'uni',0);

else % recursion - baseline_condition {unwrapped to NON-cell}, above are doing recursion (wrt to dimensions)
    

    timeBins = -conf.prior : conf.tstep : conf.after;
    % bin the spikes for this trial
    binnedSpikesOut = histcounts(spk,timeBins);
end
end

