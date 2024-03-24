session = 'T_20220602'; % session tag, 'F_20200804', 'F_20200806', 'F_20200807','R_20220304', 'R_20220316', 'R_20220318', 'T_20220526', 'T_20220602', 'T_20220603';

for curr=session
    recMonDateList = curr;
    [cum_sepa_self ,xlims, ylims]= fig1_trajsumm_mod_3d(recMonDateList, 'MN_PMv_PMd', 'self', 'active');
    cum_sepa_cross = fig1_trajsumm_mod_3d(recMonDateList, 'MN_PMv_PMd', 'cross','passive',xlims,ylims);
end



function [cum_sepa, xlims,ylims] = fig1_trajsumm_mod_3d(recMonDateList, NeuronPop, projMode, scalemode, xlims_in, ylims_in)
% ``` fig 1 traj&summary fig gen_func 
% ``` path './summ_plot/'
% input: recording_session, MN_population, projection_mode (self/cross)
% compute the cumulative separability among latend trajectories (3d ver)
%
% - recMonDateList: recording MonkeyDate List: which session (for example,
% 'T_20220602': Monkey T Date 20220602
%
% - NeuronPop: neuron population: which population (for example,
% 'MN_PMv_PMd': Mirror Neuron from all PM (PMv+PMd)
%
% - projMode: the given neural data, will be either "self" or "cross"
% projected, i.e., if 'self', the execution data will be projected into
% Execution subspaces. 
%
% - scalemode, x(y)lims_in: how the canvas is scaled. 'active' will make the canvas
% scaled by its own size. 'passive' will use the given x(y)lims. 


if ~exist('NeuronPop','var')
    NeuronPop = ['MN_', 'PMv_PMd'];   %'M1c_M1r'; PMv_PMd_M1c_M1r    PMv_PMd
end

if ~exist("projMode",'var')
    projMode = 'self';
end

if ~exist("scalemode", 'var')
    scalemode = 'active';
end


opList = {'Ex', 'Obs'};
cum_sepa = cell(1,2);

for op=1:numel(opList)
    opMode = opList{op};    
    objColorList = {'#7f58af','#64c5eb','#e84d8a','#feb326'}; % (tartget) obj color 
    
    
    behaveEpochList = {'inst','gocue','mvon','holds'}; % segments 
    subsList = {'I','G','M', 'H'}; % subspaces 
    subsIndMap = containers.Map(subsList,{"s1", "s2", "s3", "s4"});
    
    % decide the string that related to subs_time_info  -> tmpstr
    if isequal(projMode, 'self')
        tmpstr = opMode;
    elseif isequal(opMode, 'Ex')
        tmpstr = 'Obs';
    elseif isequal(opMode, 'Obs')
        tmpstr = 'Ex';
    end 
    % \nkt_data\NeuronPop_recMonDateList{1}_```tmpstr```Delay.mat', 'secs') & subs --> subspace
    % \nkt_data\NeuronPop_recMonDateList{1}_```opMode```Delay.mat',  'NKT_cond') &  --> neural_data
    
    load(strcat('./nkt_data/',NeuronPop, '_', recMonDateList,'_', tmpstr,'Delay.mat'), 'secs','subs_t');
    
    load(strcat('./nkt_data/',NeuronPop, '_', recMonDateList,'_', opMode,'Delay.mat'), 'NKT_cond');
    subsMap = containers.Map;
    for ss=subsList %  subsIndMap('I') -> subMap('I') fetch the neural_dims for each IGMH 
        subsMap(ss{:}) = subs_t{secs.(subsIndMap(ss{:})).('locind')};
    end
    
    % !note! the secs is DIFFERENT here 
    load(['./nkt_data/',NeuronPop, '_', recMonDateList,'_', opMode,'Delay.mat'], 'secs');  % use this to get neural_activity segments 
    
    % select one segm and projected onto the above subs subsMap(I|G|M|H)
    leftOffset = 0; 
    rightOffset = 99;
    secsPlotted = struct('inst', struct('h',secs.s1.locind-leftOffset, 't', secs.s1.locind+rightOffset), ...
            'gocue', struct('h',secs.s2.locind-leftOffset, 't', secs.s2.locind+rightOffset),...
            'mvon',  struct('h',secs.s3.locind-leftOffset, 't', secs.s3.locind+rightOffset),...
            'holds',  struct('h',secs.s4.locind-leftOffset, 't', secs.s4.locind+rightOffset));
     
    hfig = figure('PaperPositionMode','auto', 'Units','inches','Position',[0 0 9/2 15/2]);
    hax = gobjects(1,4*4+1); % 4x4 traj visual + 1 summ_matrix
    
    tiledlayout(6,4,"TileSpacing","tight");
    summ_mat = nan(4,4);
    
    for epochInd=1:numel(behaveEpochList)
        behaveEpoch = behaveEpochList{epochInd}; 
        ii=1; 
        for ss=subsList % inst gocue mvon holds subspace 
            % project nkt -> I,G,M,H subspaces 
            currsubs = subsMap(ss{:})'; % this sessions' I(or G,M,H) subs
            tmp1a = mean(pagemtimes(currsubs, permute(NKT_cond{1,1}(:,:,secsPlotted.(behaveEpoch).h:secsPlotted.(behaveEpoch).t), [1 3 2])),3);
            tmp2a = mean(pagemtimes(currsubs, permute(NKT_cond{1,2}(:,:,secsPlotted.(behaveEpoch).h:secsPlotted.(behaveEpoch).t), [1 3 2])),3);
            tmp3a = mean(pagemtimes(currsubs, permute(NKT_cond{1,3}(:,:,secsPlotted.(behaveEpoch).h:secsPlotted.(behaveEpoch).t), [1 3 2])),3);
            tmp4a = mean(pagemtimes(currsubs, permute(NKT_cond{1,4}(:,:,secsPlotted.(behaveEpoch).h:secsPlotted.(behaveEpoch).t), [1 3 2])),3);
            
            bias = mean(horzcat(tmp1a,tmp2a, tmp3a,tmp4a), 2); % common shifting 
            tmp1a = tmp1a - bias;
            tmp2a = tmp2a - bias;
            tmp3a = tmp3a - bias;
            tmp4a = tmp4a - bias;

            hax(1,(epochInd-1)*4+ii) = nexttile((epochInd-1)*4+ii);
            
            hold on;
            plot(tmp1a(1,:), tmp1a(2,:),'Color',objColorList{1}, 'LineWidth',2);
            plot(tmp2a(1,:), tmp2a(2,:),'Color',objColorList{2}, 'LineWidth',2);
            plot(tmp3a(1,:), tmp3a(2,:),'Color',objColorList{3}, 'LineWidth',2);
            plot(tmp4a(1,:), tmp4a(2,:),'Color',objColorList{4}, 'LineWidth',2);
            hold off;
            % add same axis 
            axis  square; 
            diff12 = tmp1a(1:3,:) - tmp2a(1:3,:);
            diff13 = tmp1a(1:3,:) - tmp3a(1:3,:);
            diff14 = tmp1a(1:3,:) - tmp4a(1:3,:);
            diff23 = tmp2a(1:3,:) - tmp3a(1:3,:);
            diff24 = tmp2a(1:3,:) - tmp4a(1:3,:);
            diff34 = tmp3a(1:3,:) - tmp4a(1:3,:);
        
            dis12 = sqrt(sum(diff12.^2,1));
            dis13 = sqrt(sum(diff13.^2,1));
            dis14 = sqrt(sum(diff14.^2,1));
            dis23 = sqrt(sum(diff23.^2,1));
            dis24 = sqrt(sum(diff24.^2,1));
            dis34 = sqrt(sum(diff34.^2,1));
        
            sumdis = dis12+dis13+dis14+dis23+dis24+dis34;
            summ_mat(epochInd, ii) = mean(sumdis);
            ii=ii+1;
        end
    end
    cum_sepa{op} = summ_mat; % cummulative separability 



    nexttile(1); % manually add scale - on the first panel 
    if (hax(16).XLim(end) -hax(16).XLim(1)) > 5
        slen = 5; % legend line_segment length 
    else
        slen = 2;
    end

    
    disp(summ_mat);
    hax(1,end) = nexttile(4*4+1,[2 4]);
    c_map  = [0,0,0; 24,24,24; 48,48,48; 96,96,96; 112,112,112; 136,136,136; 160 160 160; 176 176 176; 192 192 192; 211 211 211]./256;
    data_range = [0 1.5 2.5 4 6 8 10 20 45 55]; 
    color_mapping =  [0 1.5 2.5 4 6 8 10 20 45 55]; 
    mapped_summ_mat = interp1(data_range, color_mapping, summ_mat, 'nearest');
    imagesc(mapped_summ_mat);  % mapped_summ_mat
    colorbar;
    caxis([0 40]);
    colormap(c_map);
    for i = 1:size(summ_mat, 1)
        for j = 1:size(summ_mat, 2)
            if  isequal(summ_mat(i, j), max(summ_mat(i,:)))
                fontcolor = '#ff0000';
            else
                fontcolor = '#d3d3d3';
            end
    
            text(j, i, num2str(summ_mat(i, j), '%.1f'), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','Color',fontcolor,'FontSize',12,'FontWeight','normal');
        end
    end


    xlabel('Subspace', 'FontSize', 12);
    ylabel('Segment', 'FontSize', 12);
    axis image;
    axis on;
    grid off;
    title('  ');

    hax(end).XAxisLocation = 'bottom';
    hax(end).XTickLabel = cellfun(@(x) strcat(tmpstr(1),x), {'I', 'G', 'M','H'},'UniformOutput',false) ;
    hax(end).YTickLabel = cellfun(@(x) strcat(opMode(1),x), {'I', 'G', 'M','H'},'UniformOutput',false) ;
    set(hax(end),'Fontsize',12)

    
    
    
    % adjust the subplpot 1-16 layout, title, strings, etc..
    
    hax(1).Title.String = strcat(tmpstr(1), "I subs.");
    hax(2).Title.String = strcat(tmpstr(1), "G subs.");
    hax(3).Title.String = strcat(tmpstr(1), "M subs.");
    hax(4).Title.String = strcat(tmpstr(1), "H subs.");
    
    hax(1).YLabel.String = strcat(opMode(1), "I segm.");
    hax(1+4).YLabel.String = strcat(opMode(1), "G segm.");
    hax(1+4*2).YLabel.String = strcat(opMode(1), "M segm.");
    hax(1+4*3).YLabel.String = strcat(opMode(1), "H segm.");
    
    
    if isequal(scalemode, 'active') && isequal(opMode,'Ex') && isequal(projMode,'self')% only Ex-self use its own scale
        for i=1:16
            hax(i).XLim = deal([min(min(vertcat(hax(1:16).XLim)))  max(max(vertcat(hax(1:16).XLim)))]);
            hax(i).YLim = deal([min(min(vertcat(hax(1:16).YLim)))  max(max(vertcat(hax(1:16).YLim)))]);
        end
        xlims = [min(min(vertcat(hax(1:16).XLim)))  max(max(vertcat(hax(1:16).XLim)))];
        ylims = [min(min(vertcat(hax(1:16).YLim)))  max(max(vertcat(hax(1:16).YLim)))];
    elseif isequal(opMode,'Obs') && isequal(projMode,'self') % Obs-Self use Ex-self's scale inside the same function thus xlims
        for i=1:16
            hax(i).XLim = deal(xlims);
            hax(i).YLim = deal(ylims);
        end
    else % cross-project are using given_scale thus xlims_in
        for i=1:16
            hax(i).XLim = deal(xlims_in);
            hax(i).YLim = deal(ylims_in);
        end
    
    end
    

    
    set(hax(1:16),'xticklabel','','yticklabel','', 'FontName','Arial','FontWeight','bold', 'FontSize',12);
    set(hax(end),'FontName','Arial','FontWeight','bold', 'FontSize',12);
    
    
    exportgraphics(hfig,['./summ_plot_3d/',recMonDateList,'_', opMode,'-',projMode, '-traj&summ(', num2str(leftOffset),',', num2str(rightOffset)  ').png'],'Resolution',1000);
    close;

end




end

