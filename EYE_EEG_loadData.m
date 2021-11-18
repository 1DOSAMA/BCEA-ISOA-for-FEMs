%% EYE-EEG 信号加载，保存成标准格式
% by 杜成航 2020.10.22
% 记得先加载Tobii的SDK
%% 数据加载，同时load三个文件
[filename, pathname] = uigetfile( ...
    {'*.mat;','All Image Files';},...
    'Pick files(1)', ...
    'MultiSelect', 'on',...
    '.\DataSave');
if ischar(filename)
    files = fullfile(pathname, filename);
    load(files);
else
     for i = 1:size(filename,2)
         files = fullfile(pathname, filename{i});
         load(files);
     end
end