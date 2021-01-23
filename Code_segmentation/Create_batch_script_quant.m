function [outputArg1,outputArg2] = untitled2(inputArg1,inputArg2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
path = '/Users/denisschapiro/Dropbox (HMS)/COVID_Lung/component_images';

% Image dir
images_path = fullfile(path,'Images');
images_dir = dir([images_path '/*.tif']);
images_dir_cell = struct2cell(images_dir);
images_dir_cell_sorted = natsortfiles(images_dir_cell(1,:));

% Mask dir
mask_path = fullfile(path,'Masks');
mask_dir = dir([mask_path '/*.tiff']);
mask_dir_cell = struct2cell(mask_dir);
mask_dir_cell_sorted = natsortfiles(mask_dir_cell(1,:));


txt = {};
txt = cell(136,1);

for i=1:136
    mask_all = {};
    images_all = {};
    
    mask_all = ['/data/Masks/' mask_dir_cell_sorted{i}];
    images_all = ['/data/Images/' images_dir_cell_sorted{i}];
    
    txt{i,1} = ['python CommandSingleCellExtraction.py --masks ' mask_all ...
        ' --image ' images_all ' --output /data/ --channel_names /data/my_channels.csv'];
end
end

