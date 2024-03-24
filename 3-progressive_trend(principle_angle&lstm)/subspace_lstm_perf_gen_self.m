NeuronPop = 'nMN_PMv_PMd'; % because `n' (1+ char), below update s.t. F_20200804 - need you to change charCount by hand
Exfile_list = dir(['./nkt_data/', NeuronPop, '*Ex*.mat']);
ExRecTags = cellfun(@(x) x(8+5:17+5), {Exfile_list.name}, 'UniformOutput',false);
Obsfile_list = dir(['./nkt_data/', NeuronPop, '*Obs*.mat']);
ObsRecTags = cellfun(@(x) x(8+5:17+5), {Obsfile_list.name}, 'UniformOutput',false);
if isequal(sort(ExRecTags), sort(ObsRecTags))
    monkeyDateList = ExRecTags; % confirm Ex&Delay both exist and use 
else
    disp('data Ex/Obs not match!');
end

tstep = 50; % save time by enlarging the t_gap
iteration =10;

for mk=monkeyDateList
    MonkeyDate = mk{1}; %for example: 'F_20200807';
    load(strcat('./nkt_data/',NeuronPop, '_', MonkeyDate,'_ExDelay.mat'));
    NKT_cond_Ex = NKT_cond_lite; % use the "compressed" ver to same time, see lower where truncation_pos is divided by 10. 
    secs_Ex = secs;
    subs_t_Ex = subs_t;
    
    load(strcat('./nkt_data/',NeuronPop, '_', MonkeyDate,'_ObsDelay.mat'));
    NKT_cond_Obs = NKT_cond_lite;
    secs_Obs = secs;
    subs_t_Obs = subs_t;

    [Ex_inst_perfList_Exsubs ] = optimalSubspaceFinder(NKT_cond_Ex, subs_t_Ex, round((secs_Ex.s1.locind+0)/10), round((secs_Ex.s1.locind+100)/10),tstep,iteration);
    [Ex_gocue_perfList_Exsubs ] = optimalSubspaceFinder(NKT_cond_Ex, subs_t_Ex, round((secs_Ex.s2.locind+0)/10), round((secs_Ex.s2.locind+100)/10),tstep,iteration);
    [Ex_mvon_perfList_Exsubs ] = optimalSubspaceFinder(NKT_cond_Ex, subs_t_Ex, round((secs_Ex.s3.locind+0)/10), round((secs_Ex.s3.locind+100)/10),tstep,iteration);
    [Ex_holds_perfList_Exsubs ] = optimalSubspaceFinder(NKT_cond_Ex, subs_t_Ex, round((secs_Ex.s4.locind+0)/10), round((secs_Ex.s4.locind+100)/10),tstep,iteration);
    
    [Obs_inst_perfList_Obssubs] = optimalSubspaceFinder(NKT_cond_Obs, subs_t_Obs, round((secs_Obs.s1.locind+0)/10), round((secs_Obs.s1.locind+100)/10),tstep,iteration);
    [Obs_gocue_perfList_Obssubs] = optimalSubspaceFinder(NKT_cond_Obs, subs_t_Obs, round((secs_Obs.s2.locind+0)/10), round((secs_Obs.s2.locind+100)/10),tstep,iteration);
    [Obs_mvon_perfList_Obssubs] = optimalSubspaceFinder(NKT_cond_Obs, subs_t_Obs, round((secs_Obs.s3.locind+0)/10), round((secs_Obs.s3.locind+100)/10),tstep,iteration);
    [Obs_holds_perfList_Obssubs] = optimalSubspaceFinder(NKT_cond_Obs, subs_t_Obs, round((secs_Obs.s4.locind+0)/10), round((secs_Obs.s4.locind+100)/10),tstep,iteration);
    
    % save to the right place with formatted name
    save(['./lstm_perf/', NeuronPop, '_', MonkeyDate, '_self_',datestr(datetime,'ssmm_dd_yy'), '.mat'],...
        'Ex_inst_perfList_Exsubs','Ex_gocue_perfList_Exsubs', 'Ex_mvon_perfList_Exsubs','Ex_holds_perfList_Exsubs',...
        'Obs_inst_perfList_Obssubs','Obs_gocue_perfList_Obssubs','Obs_mvon_perfList_Obssubs','Obs_holds_perfList_Obssubs');
end




function [perfListIter, best_lstmIter, perfBestIndex] = optimalSubspaceFinder(NKT_cond, subs_t, tleft, tright,tstep, iteration)
% NKT_cond is the {1x4} cell that keeps data of each obj-task. 
% OptimalSubspace is defined as the one that ``preserves`` the best
% information such that (dim-red) data bring the best decoding performance
% given the same full-dim data. 
% % inputs:
%  - NKT_cond: 1x4 cell, each has [Neu x #Trial x Timebins] neural signal [full-dim]
%  - subs_t: subspaces, each timestamp has gen a subspace in pipline 1
%  - tleft/right: to cut neural signal @left/right 

%logic: loop thru all subspaces, for each subspace, the traj data will be
%projected (dim-red) and then the dim-red data will be used to train a
%decoder and then decoder will be tested with thre rest data. 

% classifier's job - distinguish 4 objs with given
% info-preserved-in-subspace

% steps: divide data into train/test set
% format xx set into train_nn -> use trained nn to eval perf -> get the
% best one
if nargin < 5
    tstep = 1000;
    iteration = 1;
elseif nargin < 6
    iteration = 1;
end

miniBatchSize = 10;
inputSize = 3; % feature-dimensions --> should be the subspace dim-red's dim
numHiddenUnits = 20; % 
numClasses = 4; % object(condition/task)
layers = [ ...
    sequenceInputLayer(inputSize)
    bilstmLayer(numHiddenUnits,OutputMode="last")
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];
options = trainingOptions("adam", ...
    ExecutionEnvironment="gpu", ...
    GradientThreshold=1, ...
    MaxEpochs=50, ...
    MiniBatchSize=miniBatchSize, ...
    SequenceLength="longest", ...
    Shuffle="never", ...
    Verbose=0);

tList = 1:tstep:size(subs_t,2); 
% 0325 mod - Data_Length --> subspace_length
perfListIter = zeros(iteration+1, length(tList));
perfBestIndex = zeros(iteration,1);
best_lstmIter = cell(iteration,1);
for iter=1:iteration
    ratio = 0.6; % test/all ratio
    trl1 = size(NKT_cond{1,1},2);
    labels1 = false(trl1,1); % init a flag_vector, all false indicates none of the trials is chosen
    labels1(round(trl1*ratio):trl1) = true; % choose part of the trials, shuffle
    labels1 = labels1(randperm(length(labels1)));
    
    trl2 = size(NKT_cond{1,2},2);
    labels2 = false(trl2,1); 
    labels2(round(trl2*ratio):trl2) = true;
    labels2 = labels2(randperm(length(labels2)));
    
    trl3 = size(NKT_cond{1,3},2);
    labels3 = false(trl3,1); 
    labels3(round(trl3*ratio):trl3) = true;
    labels3 = labels3(randperm(length(labels3)));
    
    trl4 = size(NKT_cond{1,4},2);
    labels4 = false(trl4,1); 
    labels4(round(trl4*ratio):trl4) = true;
    labels4 = labels4(randperm(length(labels4)));
    
    perfList = zeros(1, length(tList));
    tmpbestperf = 0;
    for ii=1:length(tList)
        if ~mod(ii,10)
            disp(strcat(num2str(ii),'timestamp'));
        end
        ti = tList(ii);
        
        tmp1 = pagemtimes(subs_t{1,ti}', permute(NKT_cond{1,1}(:,:,tleft:tright), [1 3 2]));
        tmp11 = squeeze(mat2cell(tmp1, size(tmp1,1), size(tmp1,2), ones(1,trl1)));
        tmp2 = pagemtimes(subs_t{1,ti}', permute(NKT_cond{1,2}(:,:,tleft:tright), [1 3 2]));
        tmp22 = squeeze(mat2cell(tmp2, size(tmp2,1), size(tmp2,2), ones(1,trl2)));
        tmp3 = pagemtimes(subs_t{1,ti}', permute(NKT_cond{1,3}(:,:,tleft:tright), [1 3 2]));
        tmp33 = squeeze(mat2cell(tmp3, size(tmp3,1), size(tmp3,2), ones(1,trl3)));
        tmp4 = pagemtimes(subs_t{1,ti}', permute(NKT_cond{1,4}(:,:,tleft:tright), [1 3 2]));
        tmp44 = squeeze(mat2cell(tmp4, size(tmp4,1), size(tmp4,2), ones(1,trl4)));

        trainData = [tmp11(~labels1); tmp22(~labels2); tmp33(~labels3); tmp44(~labels4)];
        trainLabel = categorical([ones(sum(~labels1),1); 2*ones(sum(~labels2),1); 3*ones(sum(~labels3),1); 4*ones(sum(~labels4),1)]);
        testData = [tmp11(labels1); tmp22(labels2); tmp33(labels3); tmp44(labels4)];
        testLabel = categorical([ones(sum(labels1),1); 2*ones(sum(labels2),1); 3*ones(sum(labels3),1); 4*ones(sum(labels4),1)]);
        
        net = trainNetwork(trainData,trainLabel,layers,options);
        YPred = classify(net,testData, MiniBatchSize=miniBatchSize, SequenceLength="longest");
        acc = sum(YPred == testLabel)./numel(testLabel);
        perfList(ii) = acc;
        if perfList(ii) > tmpbestperf
            tmpbestperf = perfList(ii); % if the perf > record, update record
            perfBestIndex(iter) = ii;
            best_lstmIter{iter} = net; % net is the nn@ current_time, update the best_lstm
        end
    end

    perfListIter(iter,:) = perfList;
end

perfListIter(end,:) = tList; % add timestamp at the very bottom row

end


