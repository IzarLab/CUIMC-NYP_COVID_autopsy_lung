function [] = CSV_analysis()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Load all csv in one table
folder_csv = ...
    '/Users/denisschapiro/Dropbox (HMS)/COVID_Lung/component_images/CSV_files';
all_csv_dir = dir([folder_csv '/*.csv']);

for i=1:size(all_csv_dir,1)
    if i==1
        table_csv = readtable(all_csv_dir(i).name);
        variable_names = table_csv.Properties.VariableNames;
        table_csv_array = table2array(table_csv);
    else
        table_csv_array = [table_csv_array;table2array(readtable(all_csv_dir(i).name))];
    end
end

% Create table
table_csv_all = array2table(table_csv_array);
table_csv_all.Properties.VariableNames=variable_names;

%% Analyze the data
% Single cell quantification comparison



end

