%% EYE-EEG �źż��أ�����ɱ�׼��ʽ
% by �ųɺ� 2020.10.22
% �ǵ��ȼ���Tobii��SDK
%% ���ݼ��أ�ͬʱload�����ļ�
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