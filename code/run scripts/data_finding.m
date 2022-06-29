% This script helps visually explore a raw recording.
% Each plot will have a unqiue color-style combination to easily identity and rule out bad electrodes.
%
%% Reset:
clear all; close all; clc;
addpath(genpath(pwd))

%% Load:
file_name = "00013547_s001_t001.edf";
[hdr, record] = edfread(file_name);

%% Signals list to ignore:
ignore_list = ["32" "31" "30" "29" "28" "27" "26" "25" "24" "23" "22" "21" "20" "PZ" "EKG1" "FZ" "CZ" "EMG" "LOC" "ROC" "PHOTIC"];
% ignore_list = [ignore_list "SP1" "SP2" "C3" "C4" "EMG" "T1" "T2" "T3" "T4" "T5" "T6" "C3P" "C4P" "P3" "P4" "O1" "O2"];
% ignore_list = [ignore_list "F7" "F8" "F3" "F4" ]
ignore_list = [ignore_list "SP1" "SP2" ];
close all;
% Plot:
plot_eeg(hdr, record, ignore_list)

%% End 

%% Subs:
function styles_and_colors = unique_styles_and_colors()
    colors = Visuals.distinguishable_colors(12);
    line_styles = {"-"; ":"; "--"}; %#ok<CLARRSTR>
    colors_cell = to_cell(colors,1);
    styles_and_colors = combinations(line_styles, colors_cell);
end

function plot_eeg(hdr, record, ignore_list)
    % Parse basic data:
    num_electrodes = size(record,1);
    labels = strings(1,num_electrodes);

    % Visual aids:
    styles_and_colors = unique_styles_and_colors();
    figure;
    progBar = Classes.PrintedProgressBar(num_electrodes, msg="Plotting...");
    styleIndex = 0;

    % Iterate:
    for i = 1 : num_electrodes
        
        progBar.step();

        % label:
        label = hdr.label{i};        
        label = StringUtils.eeg_placement(label);
        labels(i) = label;

        % Ignore unwanted signals:
        if ismember(label, ignore_list)
            continue
        end
        styleIndex = styleIndex + 1;

        % Color and style:
        try
            style = styles_and_colors{styleIndex}{1};
            color = styles_and_colors{styleIndex}{2};
        catch
            warning("Styles error. Not enough styles?");
        end

        % Plot:
        hold on
        signal = record(i,:);
        plotH = plot(signal);
        plotH.DisplayName = label;
        plotH.Color = color;
        plotH.LineStyle = style;

        % Data-Tips:
        plotH.DataTipTemplate.DataTipRows(1).Label = 'sample time';
        plotH.DataTipTemplate.DataTipRows(2).Label = 'v';
        % add:
        label_vec = repmat(label, size(signal));
        row = dataTipTextRow('label', label_vec);
        plotH.DataTipTemplate.DataTipRows(end+1) = row;

        % Labels:
        xlabel("samples")
        ylabel("v")
        title("EEG: "+string(hdr.recordID), Interpreter="none");
        grid on;
        grid minor;

        drawnow();

    end
    % End:
    progBar.clear();
    legend
    disp(hdr)
    disp("labels:");
    disp(labels)
end

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

function [str] = label_func(x,y,label)
    str = label;
end