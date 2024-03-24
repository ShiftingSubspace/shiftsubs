filelist = cellfun(@(x) strcat(x,'CC_for_plots.mat'),{'s1','s2','s3','s4'},'UniformOutput',false);
is_exlude_nmn_F = false;

if is_exlude_nmn_F
    nmn_f_startInd = 4; % exclude first three for nmn-cc's 
else
    nmn_f_startInd = 1; % default
end

% plot configs
varName = {'MN_Ex1Ex2', 'nMN_Ex1Ex2(nmn_f_startInd:end,:)',...
    'MNvnMN_Ex1','nMN_Ex1Obs1(nmn_f_startInd:end,:)','MN_Ex1Obs1',...
    'MNvnMN_Obs1','MN_Ex1Ex1','nMN_Ex1Ex1(nmn_f_startInd:end,:)'};
colorInfo = {'r','g','b','c','m','#edb120','#77ac30','k'};
lw = 2; %set ''linewidth''

figure('Units', 'normalized', 'Position', [0, 0, 1, 0.3]);
for fInd=1:numel(filelist)
    load(filelist{fInd});
    subplot(1,numel(filelist), fInd);
    for eInd=1:numel(varName)
        eval(strcat('curr=',varName{eInd},';'));
        tmp = split(varName{eInd},'('); 
        tmp = tmp{1};
        myplot(mean(curr), std(curr), struct('color',colorInfo{eInd},'lw',lw,'tag',tmp));
    end
    
    grid on;
    axis equal square;
    xlim([0,1]); xticks([0,0.5,1]);
    ylim([0,1]); yticks([0,0.5,1]);
    zlim([0,1]); zticks([0,0.5,1]);

    if isequal(fInd,1)
        xlabel('cc1'); %xticks([0,0.5,1]);
        ylabel('cc2'); 
        zlabel('cc3'); 
        legend({'MN_Ex1Ex2','','','','',...
            'nMN_Ex1Ex2', '','','','',...
            'MNvnMN_Ex1','','','','',...
            'nMN_Ex1Obs1','','','','',...
            'MN_Ex1Obs1','','','','',...
            'MNvnMN_Obs1','','','','',...
            'MN_Ex1Ex1','','','','',...
            'nMN_Ex1Ex1','','','',''},"Interpreter","none");

    else
        xticklabels([]);
        yticklabels([]);
        zticklabels([]);
    end
    
    view([30 30])
    clearvars -regexp MN*|nMN*;
end




function myplot(mu, sig, params)
if nargin < 3
    params = struct('color','r','lw',2,'tag','ccc');
end
hold on;
scatter3(mu(1),mu(2),mu(3),'Tag',params.tag,'MarkerFaceColor',params.color,'MarkerEdgeColor','none');
plot3([mu(1)-sig(1), mu(1)+sig(1)], [mu(2),mu(2)], [mu(3), mu(3)],'Color',params.color,'LineWidth',params.lw); % x-err
plot3([mu(1), mu(1)], [mu(2)-sig(2),mu(2)+sig(2)], [mu(3), mu(3)],'Color',params.color,'LineWidth',params.lw); 
plot3([mu(1), mu(1)], [mu(2),mu(2)], [mu(3)-sig(3), mu(3)+sig(3)],'Color',params.color,'LineWidth',params.lw); 


plot3([mu(1), mu(1)], [mu(2),mu(2)], [0, mu(3)],'Color','#d3d3d3','LineWidth',params.lw/2);  % marc proj line

end
% legend(findobj(gca,'Tag',params.tag), params.tag); disp(findobj('-regexp','Tag','[MN*|nMN*]'));