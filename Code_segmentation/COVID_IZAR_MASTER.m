function [] = COVID_IZAR_MASTER(folder_location)
%UNTITLED Master processing function for Covid samples from Izar lab
%   Detailed explanation goes here

%% For testing
folder_location = ...
    '/Users/denisschapiro/Dropbox (HMS)/COVID_Lung/component_images/';
% % Get image info
% info = imfinfo(image_path);

%% Process all images
filelist = dir([folder_location '*.tif']);

% Run through samples
for fileid=1:size(filelist,1)
    image_path = filelist(fileid).name;
    
    % Load images
    image_channel_1= ...
        1000*mat2gray(imread(image_path,1));
    image_channel_2= ...
        1000*mat2gray(imread(image_path,2));
    image_channel_3= ...
        1000*mat2gray(imread(image_path,3));
    image_channel_4= ...
        1000*mat2gray(imread(image_path,4));
    image_channel_5= ...
        1000*mat2gray(imread(image_path,5));
    image_channel_6= ...
        1000*mat2gray(imread(image_path,6));
    image_channel_7= ...
        1000*mat2gray(imread(image_path,7));
    image_channel_8= ...
        1000*mat2gray(imread(image_path,8));
    
    % Export image with 8 channels
    image_full = uint16(cat(3,image_channel_1,image_channel_2,image_channel_3,image_channel_4,...
        image_channel_5,image_channel_6,image_channel_7,image_channel_8));
    
    % Write image
    Ilastik_image_write(['highres_' image_path],image_full);
    
    % Create tiles
    % Random start
    random_start = [];
    random_start = round(rand*1000);
    % Create tile
    image_full_tile = [];
    image_full_tile = image_full(random_start:random_start+300,random_start:random_start+300,:);
    
    % Write image tile
    Ilastik_image_write(['tile_highres_' image_path],image_full_tile);
    
end
end

