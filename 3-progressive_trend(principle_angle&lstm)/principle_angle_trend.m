data_prefix = './nkt_data/MN_PMv_PMd_';
data_suffix ='_ExDelay.mat';

recMK = {'T_20220526','T_20220602','T_20220603'}; % replace either to generate that row % recMK = {'R_20220304','R_20220316','R_20220318'}; % recMK = {'F_20200804','F_20200806','F_20200807'};
recMonDateList = cellfun(@(x) strcat(data_prefix,x,data_suffix),recMK,'UniformOutput',false);
compMode='self'; % computation mode
saveNameTag = strcat(recMK{1}(1),data_suffix(1:3));
theta_trend(recMonDateList, compMode, saveNameTag)



function theta_trend(recMonDateList, compMode,saveNameTag)
% compute the prin_angle between reference_subspace with a sequence of
% inst_subspace

% - recMonDateList: recording MonkeyDate List: dictate which session(s) are
% plotted as a bundle. We bundle sessions from identical monkey into one,
% i.e., monkey F,R,T. 

% - compMode: computation mode: two options, 'self' or not. if used 'self',
% the reference I,G,M,H subspaces will be selected from the same task
% condition. 

% - saveNameTag: text string that names the generated figure. 

num_dataset = numel(recMonDateList);
segmentEpoch = {'s1','s2','s3','s4'}; % I,G,M,H's notations
num_bhsess = numel(segmentEpoch);
ref_subs = cell(num_dataset,num_bhsess);

if isequal(compMode,'self')
    dual_recMonDateList = recMonDateList;
else
    dual_recMonDateList =  cellfun(@(x) strrep(x, 'ExDelay', 'ObsDelay'), recMonDateList, 'UniformOutput', false);
end


for dd=1:num_dataset
    load(dual_recMonDateList{dd},"-mat","secs","subs_t");
    for ss=1:num_bhsess
        ref_subs{dd,ss} = subs_t{secs.(segmentEpoch{ss}).('locind')};
    end
    clear secs  subs_t;
end


loc_mat = nan(num_dataset, num_dataset);
lind_mat = nan(num_dataset, num_dataset);
rind_mat = nan(num_dataset, num_dataset);

for dd=1:num_dataset
    load(recMonDateList{dd},"-mat","secs");
    for ss=1:num_bhsess
        loc_mat(dd,ss) = secs.(segmentEpoch{ss}).('locind');
        lind_mat(dd,ss) = secs.(segmentEpoch{ss}).('lind');
        rind_mat(dd,ss) = secs.(segmentEpoch{ss}).('rind');
    end
    clear secs;
end

left_len = min(loc_mat - lind_mat);
right_len = min(rind_mat-loc_mat);
IGMH_pos = loc_mat(1,:);

% below will be repeated for each reference

for ss=1:num_bhsess
    % allocation memory for each reference (I,G,M,H)~(s1,s2,s3,s4)
    s1_theta_mat = nan(num_dataset, right_len(1) +left_len(1)+1);
    s2_theta_mat = nan(num_dataset, right_len(2) +left_len(2)+1);
    s3_theta_mat = nan(num_dataset, right_len(3) +left_len(3)+1);
    s4_theta_mat = nan(num_dataset, right_len(4) +left_len(4)+1);
    
    
    % function input assume chosen reference - repeat for four reference
    
    for dd=1:num_dataset
        load(recMonDateList{dd},"-mat","subs_t","secs");
        ref_sub = ref_subs{dd,ss};
        i = 1;
        for tt=loc_mat(dd,1)-left_len(1):loc_mat(dd,1)+right_len(1)  % s1
            thetas = tmp_func_theta(ref_sub, subs_t{tt});
            s1_theta_mat(dd,i) = thetas(1);
            i = i+1;
        end
    
        i = 1;
        for tt=loc_mat(dd,2)-left_len(2):loc_mat(dd,2)+right_len(2)  % s2
            thetas = tmp_func_theta(ref_sub, subs_t{tt});
            s2_theta_mat(dd,i) = thetas(1);
            i = i+1;
        end
    
        i = 1;
        for tt=loc_mat(dd,3)-left_len(3):loc_mat(dd,3)+right_len(3)  % s3
            thetas = tmp_func_theta(ref_sub, subs_t{tt});
            s3_theta_mat(dd,i) = thetas(1);
            i = i+1;
        end
    
        i = 1;
        for tt=loc_mat(dd,4)-left_len(4):loc_mat(dd,4)+right_len(4)  % s4
            thetas = tmp_func_theta(ref_sub, subs_t{tt});
            s4_theta_mat(dd,i) = thetas(1);
            i = i+1;
        end
    end
    
    
    
    hfig = figure('PaperPositionMode','auto', 'Units','inches','Position',[0 0 3.5 2.5]);
    hax = gobjects(1,13);
    lwidth=1;
    clist = {'red','green','blue'}; % can add your colors 
    tiledlayout(2,4,"TileSpacing","tight");    
    hold on;
    
    for mk_ind=1:num_dataset % each session, plot by four stages, zoom-in to see the gap
        plot(IGMH_pos(1)-left_len(1)+1:IGMH_pos(1)+right_len(1)-1, s1_theta_mat(mk_ind,2:end-1), 'Color',clist{1+mod(mk_ind,len(clist))},'LineWidth',lwidth);
        plot(IGMH_pos(2)-left_len(2)+1:IGMH_pos(2)+right_len(2)-1, s2_theta_mat(mk_ind,2:end-1), 'Color',clist{1+mod(mk_ind,len(clist))},'LineWidth',lwidth);
        plot(IGMH_pos(3)-left_len(3)+1:IGMH_pos(3)+right_len(3)-1, s3_theta_mat(mk_ind,2:end-1), 'Color',clist{1+mod(mk_ind,len(clist))},'LineWidth',lwidth);
        plot(IGMH_pos(4)-left_len(4)+1:IGMH_pos(4)+right_len(4)-1, s4_theta_mat(mk_ind,2:end-1), 'Color',clist{1+mod(mk_ind,len(clist))},'LineWidth',lwidth);
    end

    
    xline(IGMH_pos,'LineWidth',1.5);
    xticks(IGMH_pos);
    xticklabels(["I", "G","M","H"]);
    ylim([0 90]);
    yticks([0 45 90]);
    yticklabels(["0\circ", "45\circ","90\circ"])
    if ~exist('trend_plot_theta','dir')
        mkdir('trend_plot_theta');
    end
    exportgraphics(hfig,['./trend_plot_theta/', saveNameTag,'-', segmentEpoch{ss},'-',compMode,'.png'],'Resolution',600);

end


end


function thetas = tmp_func_theta(ref_sub, used_sub)
% assist func that computes the principle angle between reference subspace
% and used subspace.
[~,s,~] = svd(ref_sub'*used_sub);
thetas = real(acos(diag(s)) *180/pi);
thetas = min(thetas, 180-thetas);
end
