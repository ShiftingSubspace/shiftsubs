% visualize the CI/CD components in ``ONE`` chosen static subspace, which
% facilities the notion of "dynamical" subspace since the CD presents
% complex (tangling) patterns.

mklist = {'R_20220316'}; % 'F_20200804', 'F_20200806', 'F_20200807', 'R_20220304', , 'R_20220318', 'T_20220526', 'T_20220602', 'T_20220603'
opmode = 'ExDelay';



mk = mklist{1};
filepath = strcat('.\nkt_data\MN_PMv_PMd_',mk,'_',opmode,'.mat');
load(filepath, 'NKT_cond','secs');

tensor1 = permute(NKT_cond{1},[3 1 2]);
tensor2 = permute(NKT_cond{2},[3 1 2]);
tensor3 = permute(NKT_cond{3},[3 1 2]);
tensor4 = permute(NKT_cond{4},[3 1 2]);
tensor_all = cat(3, tensor1, tensor2, tensor3, tensor4);

numNeurons = size(NKT_cond{1},1);
numTimepoints = size(NKT_cond{1}, 3);
numTrials = cellfun(@(x) size(x,2), NKT_cond, 'UniformOutput',true);

tmp1 = nan(numTimepoints*numTrials(1),numNeurons);
for trl=1:numTrials(1)
    tmp1(1+(trl-1)*numTimepoints:trl*numTimepoints,:) = tensor1(:,:,trl);
end

tmp2 = nan(numTimepoints*numTrials(2),numNeurons);
for trl=1:numTrials(2)
    tmp2(1+(trl-1)*numTimepoints:trl*numTimepoints,:) = tensor2(:,:,trl);
end

tmp3 = nan(numTimepoints*numTrials(3),numNeurons);
for trl=1:numTrials(3)
    tmp3(1+(trl-1)*numTimepoints:trl*numTimepoints,:) = tensor3(:,:,trl);
end

tmp4 = nan(numTimepoints*numTrials(4),numNeurons);
for trl=1:numTrials(4)
    tmp4(1+(trl-1)*numTimepoints:trl*numTimepoints,:) = tensor4(:,:,trl);
end



% stack all trials into a "omni" matrix to determine the Gerneral Time
% (GT) - subspace, i.e., CI "static" subspace
tmp_all = vertcat(tmp1,tmp2,tmp3,tmp4);
[a,~,~,~,~,~] = pca(tmp_all); % copy from ../1-nexFile2nkt_data/nexmat2nktcell i.e., a=subspace, d=cumsum(d)
usedSubs = a(:,1:2);


tensor_all_mu = mean(tensor_all, 3); % avg cross all trials=sum(trl_1,2,3,4)
tensor1_mu_trl = mean(tensor1,3); % avg cross trial:obj1
tensor2_mu_trl = mean(tensor2,3); % mu cross trial:obj2
tensor3_mu_trl = mean(tensor3,3); 
tensor4_mu_trl = mean(tensor4,3);


% project R(N)->R(2)
traj_all = tensor_all_mu * usedSubs; % project mean(all_trial) -> CI 
traj1 = tensor1_mu_trl * usedSubs; % project get "orig traj = CI+CD" :obj1
traj2 = tensor2_mu_trl * usedSubs;
traj3 = tensor3_mu_trl * usedSubs;
traj4 = tensor4_mu_trl * usedSubs;

tmp_stack = vertcat(traj1,traj2,traj3,traj4,traj_all);
axis_uplims = max(tmp_stack,[],1); 
axis_lowlims = min(tmp_stack,[],1); 

% plotting 
hfig = figure(); hold on; % plot the CI (black), Sphere, Button, Coax, Perp 
multi_segment_plot(traj1, secs, struct('color','#7f58af','linewidth',2));
multi_segment_plot(traj2, secs, struct('color','#64c5eb','linewidth',2));
multi_segment_plot(traj3, secs, struct('color','#e84d8a','linewidth',2));
multi_segment_plot(traj4, secs, struct('color','#feb326','linewidth',2));
multi_segment_plot(traj_all, secs, struct('color','black','linewidth',3));
 
axis equal square % layout adjustment
xlim(1.1*[axis_lowlims(1), axis_uplims(1)]);
ylim(1.1*[axis_lowlims(2), axis_uplims(2)]);
xticks([ceil(axis_lowlims(1)), floor(axis_uplims(1))]);
yticks([ceil(axis_lowlims(2)), floor(axis_uplims(2))])


exportgraphics(hfig,['./addon_adam_CICD_traj/', mk,'-',opmode, '-CI.png'],'Resolution',1000);


% per OneNote demo: |CD|^2 << {CI}^2 -> thus need to find YAAW subspace
tmp1 = nan(numTimepoints*numTrials(1),numNeurons);
for trl=1:numTrials(1)
    tmp1(1+(trl-1)*numTimepoints:trl*numTimepoints,:) = tensor1(:,:,trl) - tensor_all_mu; % remove CI from orig
end


tmp2 = nan(numTimepoints*numTrials(2),numNeurons);
for trl=1:numTrials(2)
    tmp2(1+(trl-1)*numTimepoints:trl*numTimepoints,:) = tensor2(:,:,trl)- tensor_all_mu; % remove CI from orig;
end

tmp3 = nan(numTimepoints*numTrials(3),numNeurons);
for trl=1:numTrials(3)
    tmp3(1+(trl-1)*numTimepoints:trl*numTimepoints,:) = tensor3(:,:,trl)- tensor_all_mu; % remove CI from orig;
end

tmp4 = nan(numTimepoints*numTrials(4),numNeurons);
for trl=1:numTrials(4)
    tmp4(1+(trl-1)*numTimepoints:trl*numTimepoints,:) = tensor4(:,:,trl)- tensor_all_mu; % remove CI from orig;
end



% stack all trials determine the CD subspace, similar to above CI
tmp_all = vertcat(tmp1,tmp2,tmp3,tmp4);
[a,b,c,e,d,f] = pca(tmp_all); 
usedSubs_cd = a(:,1:2);

%     traj1cd = (tensor1_mu_trl -tensor_all_mu) * usedSubs_cd;
%     traj2cd = (tensor2_mu_trl -tensor_all_mu) * usedSubs_cd;
%     traj3cd = (tensor3_mu_trl -tensor_all_mu)* usedSubs_cd;
%     traj4cd = (tensor4_mu_trl -tensor_all_mu)* usedSubs_cd;

traj1cd = (tensor1_mu_trl -tensor_all_mu) * usedSubs;
traj2cd = (tensor2_mu_trl -tensor_all_mu) * usedSubs;
traj3cd = (tensor3_mu_trl -tensor_all_mu)* usedSubs;
traj4cd = (tensor4_mu_trl -tensor_all_mu)* usedSubs;


hfig = figure(); hold on;
multi_segment_plot(traj1cd, secs, struct('color','#7f58af','linewidth',2));
multi_segment_plot(traj2cd, secs, struct('color','#64c5eb','linewidth',2));
multi_segment_plot(traj3cd, secs, struct('color','#e84d8a','linewidth',2));
multi_segment_plot(traj4cd, secs, struct('color','#feb326','linewidth',2));


tmp_stack = vertcat(traj1cd,traj2cd,traj3cd,traj4cd);
axis_uplims = max(tmp_stack,[],1); 
axis_lowlims = min(tmp_stack,[],1); 

axis equal square
xlim(1.1*[axis_lowlims(1), axis_uplims(1)]);
ylim(1.1*[axis_lowlims(2), axis_uplims(2)]);
xticks([ceil(axis_lowlims(1)), floor(axis_uplims(1))]);
yticks([ceil(axis_lowlims(2)), floor(axis_uplims(2))]);

exportgraphics(hfig,['./addon_adam_CICD_traj/', mk,'-',opmode, '-CD.png'],'Resolution',1000);
close all;




function multi_segment_plot(traj, secs, params)
% assis_func to plot low-dim trajectories in a subspace and (avoid)
% artifact bump via plotting segment by segment 

    event_symbols = {'^','o','s','d'};
    if nargin < 3
        params = struct('color', 'red','linewidth', 2);
    end
    
    hold on;
    segs_name = fieldnames(secs);
    for ss=1:numel(segs_name)
        leftind = secs.(segs_name{ss}).('lind');
        eventind = secs.(segs_name{ss}).('locind');
        rightind = secs.(segs_name{ss}).('rind');
    
        plot(traj(leftind:rightind,1), traj(leftind:rightind,2), 'color', params.color, 'linewidth',params.linewidth);
        scatter(traj(eventind,1), traj(eventind,2), 50, event_symbols{ss},'MarkerFaceColor',params.color,'MarkerEdgeColor',params.color,...
            'linewidth',params.linewidth);
    end
    
    hold off;


end


