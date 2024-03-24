lstm_trend({'R_20220304','R_20220316','R_20220318'}, 'MN_PMv_PMd', 'Ex', 'self'); % {'T_20220526','T_20220602','T_20220603'}



function lstm_trend(recMonDateList, NeuronPop, opMode, projMode)
% treat each inst-subspace as a filter since each subspace favors different
% patterns, when high-dim data is projected into a low-dim subspace,
% certain info will be lost. 

% this function compute "decodable info" via evaluting how much
% object-related info can be extracted via lstm classifer, then plot the
% decoding accuracy with respect to each subspace (as time flows) 

% - recMonDateList: recording MonkeyDate List: specify which sessions will bundle
%
% - NeuronPop: neural population 'MN_PMv_PMd', stands for Mirror Neuron in
% all PM (PMv+PMd)
%
% - opMode: operation mode: how the RGM is done. 'Ex', monkey execute
%
% - projMode, projection mode: how the original high-dim neural data is
% projected. if 'self', for Ex data, it will be projected into a sequence
% of execution subspace. 

color_list_shade = {'#8b0000','#00ab00','#02386E'}; % dark R-G-B - each panel is a merging of three days

segmentEpoch = {'inst','gocue','mvon','holds'};
for segm = 1:numel(segmentEpoch)
    behaveEpoch = segmentEpoch{segm};   
    perfTitleMap = containers.Map(segmentEpoch,{"Instruction", "Go cue", "Movement", "Hold"});

    if isequal(projMode, 'self')
        tmpstr = opMode;
    elseif isequal(opMode, 'Ex')
        tmpstr = 'Obs';
    elseif isequal(opMode, 'Obs')
        tmpstr = 'Ex';
    end
    target_perflist=[opMode,'_',behaveEpoch,'_perfList_',tmpstr,'subs']; % the perfList to load 
    tmp = cell(1, numel(recMonDateList));
    for i = 1:numel(recMonDateList) % loop all session for identical monkey 
        currFile = dir(['./lstm_perf/', NeuronPop,'_',recMonDateList{i}, '_',projMode,'.mat']); % lstm_perf_self here PATH change 
        load(strcat('./lstm_perf/', currFile.name), target_perflist);
        load(['./nkt_data/',NeuronPop,'_', recMonDateList{i},'_',opMode,'Delay.mat'], 'secs');  %  neural data 
        load(['./nkt_data/',NeuronPop,'_', recMonDateList{i},'_',tmpstr,'Delay.mat'], 'subs_t'); % load wanted subspaces
        eval(strcat('perf{',num2str(i),'}=',target_perflist,';'));
        eval(['tmp{',num2str(i),'}=','secs;']); % keep each sesc info for each of recording session 
        eval(['subspace{',num2str(i),'}=','subs_t;']);
    end
    
    secsNameList = fieldnames(secs);
    truncatePerf = cell(length(recMonDateList), length(secsNameList));
    truncateLocInfo = cell(length(recMonDateList), length(secsNameList));
    perf_val = nan(numel(recMonDateList), length(secsNameList));  % keep each_session's perf_val for all IGMH's subspace
        
    % -- process raw data to stack them (start here) --- 
    for sInd = 1:length(secsNameList)
        secName = secsNameList{sInd};  % secName = 's2'; % repeat for all s's
        locinds = cellfun(@(x) x.(secName).('locind'), tmp, 'UniformOutput', true);
        lind = cellfun(@(x) x.(secName).('lind'), tmp, 'UniformOutput', true);
        rind = cellfun(@(x) x.(secName).('rind'), tmp, 'UniformOutput', true);
        minLeftGap = min(locinds - lind);
        minRightGap = min(rind - locinds);
        leftTrBound = locinds - minLeftGap; % left-boundary
        rightTrBound = locinds + minRightGap; % right-boundary
        
        for pp=1:length(recMonDateList)
            tt = perf{pp}(end,:);
            tt_ind = find(tt < rightTrBound(pp) & tt> leftTrBound(pp));
            truncatePerf{pp,sInd} = perf{pp}(1:end-1,tt_ind);
            truncateLocInfo{pp,sInd} = {leftTrBound(pp), rightTrBound(pp), locinds(pp)};
        end
    
        for i=1:numel(recMonDateList)
            perf_val(i,sInd) = mean(perf{i}(1:end-1, knnsearch(perf{i}(end,:)', tmp{i}.(['s',num2str(sInd)]).('locind'))));
        end
    end
    % extracting all sessions's perflist done 


    for sInd=1:length(secsNameList)
        getSecLen = cellfun(@(x) size(x,2), truncatePerf(:,sInd));
        if ~all(getSecLen == getSecLen(1))
            truncatePerf(:,sInd) = cellfun(@(x) imresize(x,[size(x,1), round(mean(getSecLen))]), truncatePerf(:,sInd), 'UniformOutput', false);
        end
    end
    
    
    
    % get concat_t, concat_y, which is real signal.
    concat_y = cell(1,length(secsNameList));
    concat_t = cell(1,length(secsNameList));
    xgap_shift = 75;
    tick_shiftpos = zeros(1, numel(secsNameList));
    for sInd = 1:length(secsNameList)
        y = cell2mat(truncatePerf(:,sInd)); % signal
        concat_y{sInd} = y;
        preext = mean(cellfun(@(x) x{3}-x{1}, truncateLocInfo(:,sInd), 'UniformOutput',true));
        postext = mean(cellfun(@(x) x{2}-x{3}, truncateLocInfo(:,sInd), 'UniformOutput',true));
        mkloc = mean(cellfun(@(x) x{3}, truncateLocInfo(:,sInd), 'UniformOutput',true));
        tt = linspace(mkloc-preext,mkloc+postext, size(y,2));
        concat_t{sInd} = tt;
        tick_shiftpos(sInd) = mkloc+xgap_shift*(sInd-1);
    end
    % -- process raw data to stack them (end here) --- 
        
    % added to report the max  value & time (optional)
    stackperf = cell2mat(concat_y);
    tmpperf = mean(stackperf);
    tmptime = cell2mat(concat_t);
    peakTime = tmptime((tmpperf == max(tmpperf)));
    peakVal = max(tmpperf);
    peakVariation = std(stackperf(:,(tmpperf == max(tmpperf))));
    initVal = tmpperf(1);
    initVariation = std(stackperf(:,1));
    if isequal(projMode, 'self') % only self-projection is valid 
        fprintf("%s: [%.2f-%.2f]->[%.2f-%.2f]   delt:(%.f, %.f, %.f, %.f)\n", behaveEpoch, initVal, initVariation, peakVal, peakVariation, peakTime - secs.s1.locind, peakTime - secs.s2.locind, peakTime - secs.s3.locind, peakTime - secs.s4.locind);
    else
        fprintf("%s: [%.2f-%.2f]->[%.2f-%.2f]   delt:(%.f, %.f, %.f, %.f)\n", behaveEpoch, initVal, initVariation, peakVal, peakVariation, 0, 0, 0, 0);
    end


    % ----------- plot section ------------- 
    
    hfig = figure('PaperPositionMode','auto', 'Units','inches','Position',[0 0 3.5 2.5]);
    hax = gobjects(1,13);
    tiledlayout(2,4,"TileSpacing","tight");
    % plot the performance curves 
    hax(1,1) = nexttile(1,[2 4]);
    sessTrlNum=size(concat_y{1},1)/numel(recMonDateList);
    for i=1:numel(recMonDateList)
        for jj=1:length(secsNameList)
            shadedErrorBar3(xgap_shift*(jj-1)+concat_t{jj}, concat_y{:,jj}(1+sessTrlNum*(i-1):sessTrlNum*i,:), color_list_shade{i},...
                1,':', 0.3, false);
        end
    end
    
    
    for sInd = 1:length(secsNameList)
        shadedErrorBar3(xgap_shift*(sInd-1)+ concat_t{sInd},concat_y{sInd},'#000000', 2, '-', 0.3, false); % plot perf-list of Sum-all
    end

    % aesthetic things 
    xline(tick_shiftpos); % added on 06072023, vline
    ylim([0 1]);
    yticks([0 0.5 1]);
    xticks(tick_shiftpos);
    xticklabels(["I", "G","M","H"]);
    hold off;
    title(strcat(perfTitleMap(behaveEpoch), ' '), 'HorizontalAlignment','left', 'FontSize',14, 'FontWeight','bold','Position', [hax(1).XLim(1) hax(1).YLim(2)]);
    % ylabel("Classification Accuracy");
    box off;
    set(hax(1), 'FontSize', 12); % all font for hax(1) - perf curve

    % purple & orange bars added by marc - 10/13/2023
    rectangle('Position',[0, 0, 500, 1/100],'FaceColor',[0.5 0 .5], 'EdgeColor',[0.5 0 .5],'LineWidth',1.5)  % purple->inst(500ms)
    rectangle('Position',[tick_shiftpos(segm), 0.99, 100, 1/100],'FaceColor',[1, 0.5 0], 'EdgeColor',[1, 0.5 0],'LineWidth',1.5);  % orange-> neu_decode_input(100ms)
    
    % port out gen'ed fig
    if ~exist("trend_plot",'dir')
        mkdir("trend_plot");
    end
    exportgraphics(hfig,['./trend_plot/', recMonDateList{1}(1),'-', behaveEpoch,'-',opMode,'-',projMode, '.png'],'Resolution',1000);


end







end