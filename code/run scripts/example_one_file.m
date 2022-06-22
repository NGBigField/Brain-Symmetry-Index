%% Reset:
clear all; close all; clc;
addpath(genpath(pwd))

%% Pretty prints:
colors = distinguishable_colors(11);
line_styles = {"-"; ":"; "--"}; %#ok<CLARRSTR> 
colors_cell = to_cell(colors,1);
styles_and_colors = combinations(line_styles, colors_cell);

%%
fullpath = "00003097_s002_t001.edf";
[hdr, record] = edfread(fullpath);
num_electrodes = size(record,1);
%%
labels = strings(1,num_electrodes);
figure;
for i = 1 : num_electrodes
    % label:
    label = hdr.label{i};
    label = label(4:end-2);
    labels(i) = label;

    % Color and style:
    style = styles_and_colors{i}{1};
    color = styles_and_colors{i}{2};

    % Plot:
    hold on
    plotH = plot(record(i,:));
    plotH.DisplayName = label;
    plotH.Color = color;
    plotH.LineStyle = style;
    
end
legend
disp(hdr)
disp("labels:");
disp(labels)
%%


%% Subs:
function [res] = combinations(vec1, vec2)
   arguments
       vec1 (:,1) cell
       vec2 (:,1) cell
   end

   L1 = length(vec1);
   L2 = length(vec2);
   cell_mat = cell(L1,L2);
   for i1 = 1 : L1
       for i2 = 1 : L2
           cell_mat{i1,i2} = {vec1{i1}, vec2{i2}};
       end
   end

   res = cell_mat(:);
end

function res = to_cell( mat, dim)
    arguments
        mat
        dim (1,1) int32 = 1
    end

    L = size(mat,dim);
    res = cell(L,1);
    for i = 1 : L 
        vec = mat(i,:);
        res{i} = vec;
    end
end