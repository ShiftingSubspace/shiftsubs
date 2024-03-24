% batch pipeline [1] - convert *.nex into NKT_cond{1,4} data and instantneous
% subspace for each timestamp based NeuronType 


NeuronType='nMN'; % nMN or MN to extract info

cortexList = {'PMv_PMd'}; % 0911 - for journal PMMN  {'PMv', 'PMd', 'PMd_PMv', 'M1r', 'M1c', 'M1r_M1c'};
addpath(genpath('./FML_matlab/'));
fileList = dir('./nexdata/*.mat');
firingAlign = 'inston'; % extraction configs: align_index, smooth_gauss
conf = struct('prior', 1.0, 'after', 4.5, 'tstep', 1/1000); % MUST be long enough to cover all(or max) length trial
width = round(50/(conf.tstep * 1000));
Nends = round(0.2/conf.tstep);
gaussian = normpdf(-Nends:Nends,0,width)/conf.tstep;

epochStatDesp = struct('e1', struct('left', 'instoff', 'right','gocue'),...
    'e2', struct('left', 'gocue', 'right', 'mvon'),...
    'e3', struct('left', 'holds', 'right', 'holde'));


for f_pt=1:length(fileList)
    nexfile = fileList(f_pt).name;
    nexFilePath = strcat('./nexdata/', nexfile); % deal with history trash - variable name not always identical 
    variableInfo = who('-file', strcat('./nexdata/',nexfile));
    if ismember('NexData', variableInfo)
        load(strcat('./nexdata/',nexfile),'NexData');
    elseif ismember('tmptmp', variableInfo)
        load(strcat('./nexdata/',nexfile),'tmptmp');
        NexData = tmptmp;
    end
    
    spikeClassifyFile = regexprep(strcat('./nexdata/',nexfile), '.mat', '_MNscreenCompiled.xls'); % MN compile from Marc 
    for op_mode = ["ExDelay", "ObsDelay"]
        [patterns, markerDict] = patternGen(protocolDecider(nexfile),op_mode);
        tMarkers = fieldnames(markerDict);
        tMarkers = tMarkers([1:5,7,8],1); % merge "5-6-7"
        for m=1:length(tMarkers)
            mkIdDict.(tMarkers{m}) = patterns{1,1}.find_marker(markerDict.(tMarkers{m}), 'reg');
        end
        sfea_options.AlignAtIndex = patterns{1,1}.find_marker(markerDict.(firingAlign)); % sfea=SpikesFindExtractAlign, 3rd index to align (``245``)
        sfea_options.AlignEnable = true; 
        [~, MT, ~,~] = SpikesFindExtractAlign(NexData.neurons{1}.timestamps,...
            NexData.events(:,2), NexData.events(:,1), patterns, sfea_options);
        % 0228 mod - secs is based on stat computation 
        epochlenstat = epochLenReport(MT, epochStatDesp, mkIdDict);
        secs.('s1') = struct('loc', 'inston', 'l', 0.5, 'r',0.5);
        secs.('s2') = struct('loc', 'gocue', 'l', epochlenstat(1), 'r', epochlenstat(2)/2);
        secs.('s3') = struct('loc', 'mvon', 'l', epochlenstat(2)/2, 'r', 0.2);
        secs.('s4') = struct('loc', 'holds', 'l', 0.2, 'r', 0.2);
        Tbins1 = round((secs.('s1').('r') + secs.('s1').('l'))./conf.tstep);
        Tbins2 = round((secs.('s2').('r') + secs.('s2').('l'))./conf.tstep);
        Tbins3 = round((secs.('s3').('r') + secs.('s3').('l'))./conf.tstep);
        Tbins4 = round((secs.('s4').('r') + secs.('s4').('l'))./conf.tstep);
        secs.('s1').('lind') = 1;
        secs.('s1').('rind') = Tbins1;
        secs.('s2').('lind') = secs.('s1').('rind')+1;
        secs.('s2').('rind') = secs.('s1').('rind')+Tbins2;
        secs.('s3').('lind') = secs.('s2').('rind')+1;
        secs.('s3').('rind') = secs.('s2').('rind') + Tbins3;
        secs.('s4').('lind') = secs.('s3').('rind')+1;
        secs.('s4').('rind') = secs.('s3').('rind') + Tbins4;
    
        secs.('s1').('locind') = secs.('s1').('lind') + round((secs.('s1').('l') + 0)./conf.tstep);
        secs.('s2').('locind') = secs.('s2').('lind') + round((secs.('s2').('l') + 0)./conf.tstep);
        secs.('s3').('locind') = secs.('s3').('lind') + round((secs.('s3').('l') + 0)./conf.tstep);
        secs.('s4').('locind') = secs.('s4').('lind') + round((secs.('s4').('l') + 0)./conf.tstep);
        
        for ctxInd =1:length(cortexList) % generate for different populations
            currentCortex = cortexList{ctxInd};
            if isequal(NeuronType, 'MN')
                valneuInd = neuIndSel(spikeClassifyFile, currentCortex);
            elseif isequal(NeuronType, 'nMN')
                [~,valneuInd] = neuIndSel(spikeClassifyFile, currentCortex);
            end
            NKT_cond = cell(1,length(patterns)); %[N]euron, [K]trials, [T]imebins
            for cond=1:length(patterns)
                % 0228 mod 
                tmpTime = cell2mat(MT{cond});
                % get stat based on epochDesp 
                % use epoch Desp to gen secs struct & assign Tbins1+Tbins2+Tbin3 
                NKT_mat = zeros(length(valneuInd), length(MT{cond}), Tbins1+Tbins2+Tbins3+Tbins4); % Tbins change 
                for n=1:length(valneuInd)
                    disp(strcat(" N: ", num2str(n), '/', num2str(length(valneuInd))));
                    nInd = valneuInd(n);
                    spikes_ts =  SpikesFindExtractAlign(NexData.neurons{nInd}.timestamps,...
                            NexData.events(:,2), NexData.events(:,1), patterns{1,cond}, sfea_options);
                    spikes_ts = spikes_ts{1,1};
                    binnedSpikes = dpcabin_formatSpikes(spikes_ts, conf); 
                    rate = cellfun(@(x) sqrt(conv(x, gaussian, 'same')), binnedSpikes, 'UniformOutput', false);
                    for trl=1:length(rate)
                        currRate = rate{trl,1};
                        currMarkerTime = MT{1,cond}{trl,1};
                        for ss=1:length(fieldnames(secs))
                            thissec = strcat('s', num2str(ss));
                            anc_t = currMarkerTime(mkIdDict.(secs.(thissec).('loc'))) - currMarkerTime(sfea_options.AlignAtIndex);
                            left_t = anc_t - secs.(thissec).('l');
                            right_t = anc_t + secs.(thissec).('r');
            
                            leftIndex = round((conf.prior + left_t)./conf.tstep);
                            rightIndex = leftIndex + length(secs.(thissec).('lind'):secs.(thissec).('rind')) - 1; % change to this passive length 
                            NKT_mat(n,trl,secs.(thissec).('lind'):secs.(thissec).('rind')) = currRate(leftIndex:rightIndex);
                        end
                    end
                end
                NKT_cond{cond} = NKT_mat;
            end
            
            subs_t = cell(1,secs.s4.rind); 
            subs_cumvar = cell(1,secs.s4.rind);
            for tt=1:secs.s4.rind 
                tmp1 = mean(squeeze(NKT_cond{1}(:,:,tt)),2); % use avg(trials) [adam] or not - debug idea! - geometry theorem proof by Marc 2dots->3dots->4dots
                tmp2 = mean(squeeze(NKT_cond{2}(:,:,tt)),2);
                tmp3 = mean(squeeze(NKT_cond{3}(:,:,tt)),2);
                tmp4 = mean(squeeze(NKT_cond{4}(:,:,tt)),2);
                tmp = [tmp1, tmp2, tmp3, tmp4]'; % here to use four obj together to define a subspace
                [a,b,c,e,d,f] = pca(tmp); % a:=(subspace coordinate coeff) c(var)
                subs_t{tt} = a;
                subs_cumvar{tt} = cumsum(d); % 0911-mod-add cum_var
            end % 
        
            % store the NKT_cond, subs_t variables with self-explained fileName 
            if ~exist('nkt_data', 'dir')
                mkdir("./nkt_data");
            end
            fprintf(strcat("./nkt_data/",NeuronType,'_', currentCortex,'_', nexfile(1:10),'_', op_mode,'.mat',' is to generated....\n'));
            save(strcat("./nkt_data/",NeuronType,'_', currentCortex,'_', nexfile(1:10),'_', op_mode,'.mat'),'subs_t','subs_cumvar','NKT_cond','secs');
        end
    end
end
run("./caller_nkt_Tsqueezer.m");

% ----- memos ----
% Output: NeuronType_Cortex_Monkey_Date_Operation.mat 
%          MN          PMv    F     20200804 Ex
% nex -> NKD data formatter job completed 
% followed with new2_main_aim12_combined.m (from pipline2.m) 