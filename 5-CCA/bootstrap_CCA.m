% script that generates booostrapping CCA results, which is later used in
% ``./python_cca_visual/caller.ipynb`` % ``divergent_sim_allinone.mergefigs("GENERATED_MAT_FILE")``
% INPUTS: 
%   - : specify which population 
%   - iter: # iteration to bootstrap (500)
%   - truncates: specify which truncation used to compute 
%   - rec1/2 - which recording sessions are used to align 


NeuronPop = 'MN_PMv_PMd'; % Neuron Population
iter = 500; % boot-iter 
numTrialUsed = 20; % in each iter, how many trials are selected (randomly, with replacement)
truncateLeftOffset = 0; % head of truncation_position
truncateRightOffset = 99; % tail of truncation_position


reclist      = {'T_20220602'}; % {'F_20200804', 'F_20200804', 'F_20200806', 'R_20220316', 'T_20220526', 'T_20220526', 'T_20220603'};
reclist_diff = {'T_20220603'};  % {'F_20200806', 'F_20200807', 'F_20200807', 'R_20220318', 'T_20220602', 'T_20220602', 'T_20220603'};

for rr=1:numel(reclist)
    rec1=reclist{rr};
    rec2 =reclist_diff{rr};

sessSymbol = 's1'; % instruction 
% Ex1 vs Obs1
CCA_index_Ex1vObs1_s1 = cca_bootstrap(NeuronPop, sessSymbol, iter, [rec1,',Ex'], [rec1,',Obs'], numTrialUsed, truncateLeftOffset, truncateRightOffset);
% Ex1 vs Ex2
CCA_index_Ex1vEx2_s1 = cca_bootstrap(NeuronPop, sessSymbol, iter, [rec1,',Ex'], [rec2,',Ex'], numTrialUsed, truncateLeftOffset, truncateRightOffset);




sessSymbol = 's2'; % gocue
% Ex1 vs Obs1
CCA_index_Ex1vObs1_s2 = cca_bootstrap(NeuronPop, sessSymbol, iter, [rec1,',Ex'], [rec1,',Obs'], numTrialUsed, truncateLeftOffset, truncateRightOffset);
% Ex1 vs Ex2
CCA_index_Ex1vEx2_s2 = cca_bootstrap(NeuronPop, sessSymbol, iter, [rec1,',Ex'], [rec2,',Ex'], numTrialUsed, truncateLeftOffset, truncateRightOffset);



sessSymbol = 's3'; % mvmt
% Ex1 vs Obs1
CCA_index_Ex1vObs1_s3 = cca_bootstrap(NeuronPop, sessSymbol, iter, [rec1,',Ex'], [rec1,',Obs'], numTrialUsed, truncateLeftOffset, truncateRightOffset);
% Ex1 vs Ex2
CCA_index_Ex1vEx2_s3 = cca_bootstrap(NeuronPop, sessSymbol, iter, [rec1,',Ex'], [rec2,',Ex'], numTrialUsed, truncateLeftOffset, truncateRightOffset);


sessSymbol = 's4'; % hold
% Ex1 vs Obs1
CCA_index_Ex1vObs1_s4 = cca_bootstrap(NeuronPop, sessSymbol, iter, [rec1,',Ex'], [rec1,',Obs'], numTrialUsed, truncateLeftOffset, truncateRightOffset);
% Ex1 vs Ex2
CCA_index_Ex1vEx2_s4 = cca_bootstrap(NeuronPop, sessSymbol, iter, [rec1,',Ex'], [rec2,',Ex'], numTrialUsed, truncateLeftOffset, truncateRightOffset);

save(strcat("./python_cca_visual/", NeuronPop, '-',rec1,'-',rec2,'-', num2str(iter),'.mat'), 'CCA_index_*');

end



function [CCA_index] = cca_bootstrap(NeuronPop, sessSymbol, iter, sampleA, sampleB, numTrialUsed, truncateLeftOffset, truncateRightOffset)
% boostrap of CCA, use MonkeyDateA,B to check similarity between two trajs from different sessions 

if ~exist("numTrialUsed","var")
    numTrialUsed = 20;
end

if ~exist("truncateLeftOffset","var")
    truncateLeftOffset = 20;
end

if ~exist("truncateRightOffset","var")
    truncateRightOffset = 199;
end
% parsing file name
sess2Name = containers.Map({'s1','s2','s3','s4'}, {'_inst_perfList_', '_gocue_perfList_', '_mvon_perfList_', '_holds_perfList_'}); % sess2Name('s1') dict s->perfList char
CCA_index = nan(3,iter); % each column -> cca_index for that iteration iter:#bootstrap

infoA = strsplit(sampleA,','); %```F_20200806,Ex``` format example
MonkeyDateA = infoA{1}; 
OpmodeA = infoA{2};

infoB = strsplit(sampleB,',');
MonkeyDateB = infoB{1};
OpmodeB = infoB{2};

% find the ``optimalSubsapce`` based on perf_lists for each session 
% then will project high-dim SinConcatMat onto the related subspace (data
% can change as randomly choosing trials, but subspaces stays the same)
load(['./nkt_data/', NeuronPop,'_', MonkeyDateA, '_', OpmodeA, 'Delay.mat'],'subs_t');
optSubA = subs_t{optimalSubsIndFinder(['./lstm_perf/', NeuronPop, '_', MonkeyDateA, '_self.mat'], [OpmodeA, sess2Name(sessSymbol),OpmodeA,'subs'])};
load(['./nkt_data/', NeuronPop,'_', MonkeyDateB, '_', OpmodeB, 'Delay.mat'],'subs_t');
optSubB = subs_t{optimalSubsIndFinder(['./lstm_perf/', NeuronPop, '_', MonkeyDateB, '_self.mat'], [OpmodeB, sess2Name(sessSymbol),OpmodeB,'subs'])};
% confirmed *self* is used


for ii=1:iter
    SinTrlConcateMatrixA=sinTrlConcatMatrix_CCA_CD(['./nkt_data/', NeuronPop,'_', MonkeyDateA, '_', OpmodeA, 'Delay.mat'], sessSymbol, numTrialUsed, truncateLeftOffset,truncateRightOffset);
    
    SinTrlConcateMatrixB=sinTrlConcatMatrix_CCA_CD(['./nkt_data/', NeuronPop,'_', MonkeyDateB, '_', OpmodeB, 'Delay.mat'], sessSymbol, numTrialUsed, truncateLeftOffset,truncateRightOffset);
    
    LatentSinTrlConTrajA = optSubA'*SinTrlConcateMatrixA;
    LatentSinTrlConTrajB = optSubB'*SinTrlConcateMatrixB;
    [~,~,~,S,~]=latentTrajAligner_CCA(LatentSinTrlConTrajA,LatentSinTrlConTrajB); % temp check 
    CCA_index(:,ii) = diag(S);
end

end



function optimalSubsIndex=optimalSubsIndFinder(path2perf, perflistName)
    % find the index of the optimal subspace from the given perfList 
    load(path2perf,perflistName);
    eval(['currPerfList =', perflistName,';']);
    y = mean(currPerfList(1:end-1,:),1);
    optimalSubsIndex = currPerfList(end, (y==max(y)));
end


function SinTrlConcateMatrix=sinTrlConcatMatrix_CCA_CD(path2file, interestSecSYm, numTrialUsed, truncateLeftOffset,truncateRightOffset)
    % single trial concatnated matrix generator. 
    % use `path2file` to truncate the given `interestSecSYm` ('s1' to
    % truncate around instrction, 's2' around gocue, etc) for
    % `numTrialUsed` (for example, 20) times with specified
    % truncation_position in truncateLeft(Right)Offset

    % this guarantees that the chosen sampled concat-matrix is always of
    % the same size, which is required by CCA (see Juan Gallego's paper for
    % details).

    load(path2file,'NKT_cond_CD'); % N-k-T data to be truncated, concatenated
    load(path2file,'secs'); % location specifier
    
    truncatePosition = secs.(interestSecSYm).locind;
    truncateLen = truncateRightOffset+truncateLeftOffset+1;
    
    neuNum = size(NKT_cond_CD{1},1);
    objNum = size(NKT_cond_CD,2);
    objTrialNum = cellfun(@(x) size(x,2), NKT_cond_CD,'UniformOutput',true);
    if min(objTrialNum) < numTrialUsed % exception dealing with #trials
        disp("asked num2TrialUsed is more than the recorded trials, use the acutal min of trials");
        numTrialUsed = min(objTrialNum);
    end
    
    truncateSignal = cell(1,objNum);
    for oo=1:objNum
        tmp = nan(neuNum, truncateLen*numTrialUsed);
        thisObjTrail = objTrialNum(oo);
        trialSelected = randperm(thisObjTrail, numTrialUsed);
        for trl=1:numTrialUsed
            thistrl = trialSelected(trl);
            tmp(:,truncateLen*(trl-1)+1:truncateLen*trl) = NKT_cond_CD{1,oo}(:,thistrl,truncatePosition-truncateLeftOffset:truncatePosition+truncateRightOffset);
        end
        truncateSignal{oo} = tmp;
    end
    SinTrlConcateMatrix = cell2mat(truncateSignal);
end


function [Ma,Mb,U,S,V]=latentTrajAligner_CCA(La,Lb)
    % make sure the La/Lb time, dimension are identical 
    % gallergo papers uses 10/8 for different monkeys, but the
    % time=timestamp*trial*target is always identical across
    % sessions(day1,day2)
    
    % make sure La each ``row`` stands for a dimension La/Lb should be (dim x T)
    latentDim = size(La,1);
    [qa,ra] = qr(La');
    [qb,rb] = qr(Lb');
    
    [U,S,V] = svd(qa(:,1:latentDim)'*qb(:,1:latentDim)); % use the m-dim columns
    Ma = ra(1:latentDim,:)\U;
    Mb = rb(1:latentDim,:)\V;

end