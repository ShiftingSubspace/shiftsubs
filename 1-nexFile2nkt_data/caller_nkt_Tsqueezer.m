% NKT_Tsqueezer.m caller used to shrink the NKT's T such that each ten
% point is 1. NKT_cond_lite will be patched to the original target mat file

% Set the folder path
folder_path = './nkt_data/';

% Get a list of all files in the folder
file_list = dir(folder_path);

% Filter the list to include only MAT files
mat_file_list = {file_list(endsWith({file_list.name}, '.mat')).name};

for p=1:length(mat_file_list)
    file2deal = mat_file_list{1,p};
    fprintf(strcat(file2deal, ' patching...\n'));
    NKT_Tsqueezer(strcat(folder_path,file2deal));
end