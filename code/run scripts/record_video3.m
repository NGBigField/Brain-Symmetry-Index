%% clear
clearvars('-except','Data');
addpath(genpath(pwd));
close all; clc;

%% Constants:
% Choose if to filter:
use_pre_bandpass_filter = true;
% Define what to plot:
requestedData.placements = [];
requestedData.condition  = [];

%% Prepare:
figRecorder = Classes.FigureRecorder(); % Prepare Recorder
% Adjust plot limits:
Limits = derive_plot_limits(Data.EEGs);
Limits.time = Limits.time/2;
Limits.freq = Limits.freq/4;
Limits.freq(1) = 0;
Limits.bsi(1) = 0;
Limits.bsi(2) = 2;

%% Prepare figures and axes:
figH = figure;
figH.Position(1) = 5;
figH.Position(2) = 65;
figH.Position(3) = 2.3*figH.Position(3); % longer figure;
figH.Position(4) = 1.5*figH.Position(4); % longer figure;

% general values
Colors = colororder;
Sides = ["left", "right"];
cyclic_colors = @(i) get_cyclic_color(Colors, i);

% arange locations in grid:
m=2; n=5;
ind = @(x,y) Visuals.grid_xy_to_index(m,n,x,y);
subplotHandles = struct();
subplotHandles.time     = subplot(m,n,ind(1,1));
subplotHandles.freq     = subplot(m,n,ind(2,1));
subplotHandles.hhtLeft  = subplot(m,n,ind([1,2],2));
subplotHandles.hhtRight = subplot(m,n,ind([1,2],3));
subplotHandles.bsi      = subplot(m,n,ind([1,2],4));
subplotHandles.bAsy     = subplot(m,n,ind([1,2],5));
subplotHandles_hht_array = [subplotHandles.hhtLeft, subplotHandles.hhtRight];

% Shift subplots:
enlarge_factor = 1.95;
reduce_factor = 0.1;
Visuals.change_axis_width( subplotHandles.time, [] , enlarge_factor )
Visuals.sequencialy_change_axes_widths( [subplotHandles.freq, subplotHandles.hhtLeft, subplotHandles.hhtRight, subplotHandles.bsi, subplotHandles.bAsy], enlarge_factor );
Visuals.sequencialy_change_axes_widths( [subplotHandles.bsi, subplotHandles.bAsy] , reduce_factor );
keys = fieldnames(subplotHandles);
for iK = 1 : length(keys)
    key = string(keys(iK));
    subplotHandles.(key).Position(1) = subplotHandles.(key).Position(1) - 0.08;
end

% subplotHandles.freq.Position(4) = 0.9 * subplotHandles.freq.Position(4);

% time:
axisH = subplotHandles.time;
pretty_plot(axisH)
xlabel(axisH, "time [sec]",FontSize=15, Interpreter="none")
ylabel(axisH, "$$x(t)$$"  ,FontSize=16, Interpreter="latex")
title(axisH, "Time Analysis");
ylim(axisH, Limits.time)

% freq:
axisH = subplotHandles.freq;
pretty_plot(axisH)
xlabel(axisH, "freq [Hz]", FontSize=15, Interpreter="none")
ylabel(axisH, "$$|X^{F}(f)|$$"  ,FontSize=16, Interpreter="latex")
title(axisH, "Frequency Analysis");
ylim(axisH, Limits.freq)

% hht:
for iSide = 1 : 2
    axisH = subplotHandles_hht_array(iSide);
    side = Sides(iSide);
    pretty_plot(axisH)
    xlabel(axisH, "time [sec]", FontSize=15, Interpreter="none")
    ylabel(axisH, "freq [hz] ", FontSize=15, Interpreter="none")
    title(axisH, "HHT "+side);
    ylim(axisH, Limits.hht)
end

% bsi:
axisH = subplotHandles.bsi;
xlim( axisH, [0,1/2]);
ylim( axisH, Limits.bsi);
title(axisH, "BSI")

% bAsy:
axisH = subplotHandles.bAsy;
xlim( axisH, [0,1/2]);
ylim( axisH, Limits.bAsy);
title(axisH, "BAsy")

%% Iterate:
nData = length(Data.EEGs);
progBar = Classes.ProgressBar(nData, cancelable=true);

for iData = 1 : nData
    %% parse data
    eegData = Data.EEGs(iData);
    left  = eegData.left;
    right = eegData.right;
    time  = eegData.time;
    fs    = eegData.left.props.frequency;
    Sigs  = [left.recording, right.recording];

    %% Ignore some Data:
    isSkip = filter_data(eegData, requestedData);
    if isSkip
        progBar.step();
        continue
    end    

    %% Filter Signal:
    if use_pre_bandpass_filter
        Sigs(:,1) = DSP.band_pass_filter(Sigs(:,1), fs);
        Sigs(:,2) = DSP.band_pass_filter(Sigs(:,2), fs);
    end

    %% Title with global info:
    super_title_str = string(eegData.condition)+newline+left.placement.name+"-"+right.placement.name; 
    sgtitle(super_title_str);
    
    %% Plot sigs in time:
    axisH = subplotHandles.time;
    Visuals.clear_old_lines(axisH)
    for iSig = [1, 2]
        [sig, color, side] = deal_data(iSig, Sigs, Colors, Sides);
        % plot:
        hold( axisH, "on");
        plotH = plot( axisH, time, sig);
        line_props(plotH, side, color)
    end

    %% Plot freqs:
    axisH = subplotHandles.freq;
    Visuals.clear_old_lines(axisH)
    [L, R, f] = Algo.filteredFFTs(left, right);    
    cellItems = {{L,left},{R,right}};
    for iC = 1 : length(cellItems)
        X   = cellItems{iC}{1};
        sig = cellItems{iC}{2};
        [~, color, side] = deal_data(iC, Sigs, Colors, Sides);
        hold( axisH, "on");
        plotH = plot( axisH, f, abs(X));
        line_props(plotH, side, color)
    end

    %% Plot BSI:
    axisH = subplotHandles.bsi;
    Visuals.clear_old_lines(axisH);
    bsi = Algo.singleBSI(L,R);        
    hold(axisH, "on");
    plotH = bar( subplotHandles.bsi, 0, bsi,1, "stacked", "blue");

    %% Gather IMFs and HHT information:
    [imf_tensor, residues] = Algo.bemd(Sigs);
    num_components = size(imf_tensor, 3);    
    IMFs_to_include = get_IMFs_to_include(num_components);
    nImf = length(IMFs_to_include);
    HHTs_cell = cell(2, nImf);
    for iSide 1 : 2
        for iImf = 1 : nImf
            imf_num = IMFs_to_include(iImf);
            imf = imf_tensor(:, iSide, imf_num);
            [t, f, h] = Algo.hht(imf, fs);
        end
    end

    % also keep hht components:
    %% Plot HHT:
    for iSide = 1 : 2
        axisH = subplotHandles_hht_array(iSide);
        Visuals.clear_old_lines(axisH);
        for iImf = 1 : length(IMFs_to_include)
            imf_num = IMFs_to_include(iImf);
            color = cyclic_colors(iImf);
            imf = imf_tensor(:, iSide, imf_num);
            name = "imf-"+string(imf_num);            
            [t, f, h] = Algo.hht(imf, fs);
            hold(axisH, "on")
            scatter(axisH, t, f, 3, color, "filled", DisplayName=name)
        end
    end
    %% Capture video:
    figRecorder.capture(figH);

    %% Update Keep Track of Prgoress:
    cancel_requested = progBar.step();
    if cancel_requested
        break
    end
end

%% Close:
progBar.close();
figRecorder.finish();
Sounds.finish()
disp("End.")

%% End:


%% Subs:
function line_props(plotH, side, color) 
    arguments
        plotH (1,1) matlab.graphics.chart.primitive.Line
        side  (1,1) string
        color 
    end 
    plotH.DisplayName = side;
    plotH.Color = color;
    plotH.LineWidth = 1;
end
%%
function pretty_plot(axisH)
    legH = legend(axisH);
    legH.Interpreter = 'none';
    legH.FontSize = 12;
    % grid(axisH,"on");
end
%% 
function [Limits] = derive_plot_limits(EEGs)
    %% Init Global limits:
    Limits.time = [inf, -inf];
    Limits.freq = [inf, -inf];
    Limits.bsi  = [inf, -inf];
    Limits.hht  = [inf, -inf];
    Limits.bAsy = [inf, -inf];

    % Init crnt Limits:
    crntLims = struct;
    fieldNames = string(fieldnames(Limits));
    for iF = 1 : length(fieldNames)
        fieldName = fieldNames(iF);
        crntLims.(fieldName) = struct;
        crntLims.(fieldName).min = [];
        crntLims.(fieldName).max = [];
    end


    %% Iterate:
    nData = length(EEGs);
    for i = 1 : nData

        % get data:
        eegData = EEGs(i);
        left  = eegData.left;
        right = eegData.right;
        [L, R, ~] = Algo.filteredFFTs(left, right);      
        bsi = Algo.singleBSI(L,R);
        hht = [0, 30];
        bAsy = [-2, 2];

        % Derive current limits:
        crntLims.time.max = max([left.recording; right.recording]);
        crntLims.time.min = min([left.recording; right.recording]);
        crntLims.freq.max = max([abs(L); abs(R)]);
        crntLims.freq.min = min([abs(L); abs(R)]);
        crntLims.bsi.max  = bsi;
        crntLims.bsi.min  = bsi;
        crntLims.hht.max  = max(hht);
        crntLims.hht.min  = min(hht);
        crntLims.bAsy.max = max(bAsy);
        crntLims.bAsy.min = min(bAsy);

        % Update global limits:
        cellItems = {{"min", @min, 1}, {"max", @max, 2}};
        for iF = 1 : length(fieldNames)
            valFieldName = fieldNames(iF);

            for iC =  1 : length(cellItems)
                str    = cellItems{iC}{1};
                minmax = cellItems{iC}{2};
                index  = cellItems{iC}{3};
                Limits.(valFieldName)(index) = minmax( Limits.(valFieldName)(index) , crntLims.(valFieldName).(str) );
            end
        end % for iF

    end % for i 
end % func
%%
function  isSkip = filter_data(eegData, requestedData)
    arguments
        eegData (1,1) Classes.EEGData
        requestedData (1,1) struct
    end
    % Derive:
    pLeft  = eegData.left.placement.name;
    pRight = eegData.right.placement.name;
    % Decide:
    isSkip = false;
end
%% 
function [sig, color, side] = deal_data(i, Sigs, Colors, Sides)
    sig = Sigs(:,i);
    side = Sides(i);
    color = Colors(i, :);
end
%%
function indices = get_IMFs_to_include(N)
    arguments
        N (1,1) {mustBeInteger}
    end
    if N < 2
        error("Not enough IMFs");
    end
    first = 2;
    last = max(2, N-1);
    indices = first : last ;
end
%%
function color = get_cyclic_color(Colors, ind_in)
    N = size(Colors, 1);
    ind_used = mod(ind_in-1,N)+1;
    color = Colors(ind_used, :);
end