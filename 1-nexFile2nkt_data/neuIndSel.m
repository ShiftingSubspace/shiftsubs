function [mirrneuIndex, nonmirrIndex,MNTags,nonMNTags] = neuIndSel(unitFile, areaChosen)
% Use regex () to choose ideal target neurons 
% ref:https://regexone.com/lesson/excluding_characters
% ref:https://riptutorial.com/regex/example/5023/matching-various-numbers
if nargin < 2
    areaRegex = '(PMv$|PMd$|M1c$|M1r$)';
else
    areaRegex = strcat('(', strrep(areaChosen,'_','$|'), '$)');
end

T = readtable(unitFile);
% T = readtable("./F_20210412_RGM_ExeObs_0072_SH_sorted_spikeClassification.xls");
snr_filter = cell2mat(cellfun(@(x) ~isempty(regexpi(x, '(P$|M$|S$)', 'once')), table2cell(T(:,4)), 'UniformOutput', false));
cortex_filter = cell2mat(cellfun(@(x) ~isempty(regexpi(x, areaRegex, 'once')), table2cell(T(:,6)), 'UniformOutput', false));

pval = 0.05/6;
CI_Ex = cell2mat(cellfun(@(x) (x< pval), table2cell(T(:,8)), 'UniformOutput', false));
CI_Obs = cell2mat(cellfun(@(x) (x< pval), table2cell(T(:,9)), 'UniformOutput', false));

CD_Ex = cell2mat(cellfun(@(x) (x< pval), table2cell(T(:,14)), 'UniformOutput', false));
CD_Obs = cell2mat(cellfun(@(x) (x< pval), table2cell(T(:,15)), 'UniformOutput', false)); 

CX_Ex = cell2mat(cellfun(@(x) (x< pval), table2cell(T(:,20)), 'UniformOutput', false));
CX_Obs = cell2mat(cellfun(@(x) (x< pval), table2cell(T(:,21)), 'UniformOutput', false));

mn_filter = (CI_Ex|CD_Ex|CX_Ex) & (CI_Obs|CD_Obs|CX_Obs);
nmn_filter = (CI_Ex|CD_Ex|CX_Ex) & ~(CI_Obs|CD_Obs|CX_Obs);
% mn_filter = cell2mat(cellfun(@(x) ismember(x, [3 6 7 9 11:15 18 19 22:31 33 35:39 41 43:63]), table2cell(T(:,10)), 'UniformOutput', false));
% nmn_filter = cell2mat(cellfun(@(x) ismember(x, [2 8 10 32 34 40 42]), table2cell(T(:,10)), 'UniformOutput', false));
mirror_neuInd= snr_filter & cortex_filter & mn_filter;
nmirror_neuInd= snr_filter & cortex_filter & nmn_filter;
mirrneuIndex = find(mirror_neuInd == 1);
nonmirrIndex = find(nmirror_neuInd == 1);

% mod 05092023 - the details(#Neurons for each sub_pop)
subpopTags = table2cell(T(:,6));
MNTags = subpopTags(mirrneuIndex); % each cell includes a string "tag"
nonMNTags = subpopTags(nonmirrIndex);

end