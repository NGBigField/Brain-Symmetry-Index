%% clear
clearvars('-except','Data');
addpath(genpath(pwd));
close all; clc;

%% Define what to plot:
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
figH.Position(1) = 1;
figH.Position(2) = 1;
figH.Position(3) = 2.3*figH.Position(3); % longer figure;
figH.Position(4) = 1.4*figH.Position(4); % longer figure;

m=1; n=5;
subplotHandles = struct();
subplotHandles.time = subplot(m,n,1);
subplotHandles.freq = subplot(m,n,2);
subplotHandles.hht  = subplot(m,n,3);
subplotHandles.bsi  = subplot(m,n,4);
subplotHandles.bAsy = subplot(m,n,5);

factor = 1.8;
subplotHandles_array = struct2array(subplotHandles);
Visuals.widen_axes( subplotHandles_array, 1, factor )
Visuals.widen_axes( subplotHandles_array, 2, factor )
Visuals.widen_axes( subplotHandles_array, 3, factor )
Visuals.widen_axes( subplotHandles_array, 4, 0.1 )
Visuals.widen_axes( subplotHandles_array, 5, 0.1 )

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
axisH = subplotHandles.hht;
pretty_plot(axisH)
xlabel(axisH, "time [sec]", FontSize=15, Interpreter="none")
ylabel(axisH, "freq [hz] ", FontSize=15, Interpreter="none")
title(axisH, "HHT");
ylim(axisH, Limits.hht)

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

% general values
Colors = colororder;
Sides = ["left", "right"];

%% Iterate:
nData = length(Data.EEGs);
progBar = Classes.ProgressBar(nData, cancelable=true);

for i = 1 : nData
    %% parse data
    eegData = Data.EEGs(i);
    left  = eegData.left;
    right = eegData.right;
    time  = eegData.time;
    Sigs  = [left.recording, right.recording];

    %% Filter Data:
    isSkip = filter_data(eegData, requestedData);
    if isSkip
        progBar.step();
        continue
    end    

    %% Title with global info:
    sgtitle(string(eegData.condition));
    
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
    bsi = Algo.singleBSI(L,R);    
    Visuals.clear_old_lines(axisH);
    hold(axisH, "on");
    plotH = bar( subplotHandles.bsi, 0, bsi,1, "stacked", "blue");
    
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
    grid(axisH,"on");
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
        hht = [0, 40];
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