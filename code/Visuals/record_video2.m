%% clear
clearvars('-except','Data');
close all; clc;

%% Prepare:
figRecorder = Classes.FigureRecorder(); % Prepare Recorder
% Adjust plot limits:
Limits = derive_plot_limits(Data.EEGs);
Limits.time = Limits.time/2;
Limits.freq = Limits.freq/4;
Limits.freq(1) = 0;
Limits.bsi(1) = 0;

%% Iterate:
nData = length(Data.EEGs);
progBar = Classes.ProgressBar(nData);
for i = 1 : nData
    %% parse data
    eegData = Data.EEGs(i);
    left  = eegData.left;
    right = eegData.right;
    time  = eegData.time;

    %% Figure:
    figH = figure;
    figH.Position(3) = 2.0*figH.Position(3); % longer figure;
    figH.Position(4) = 1.4*figH.Position(4); % longer figure;
    % sgtitle(left.placement.name+"+"+right.placement.name)
    sgtitle(string(eegData.condition));
    %% Plot sigs:    
    axisH = subplot(1,3,1);
    for sig = [left, right]
        record = sig.recording;
        hold(axisH, "on");
        plotH = plot(time, record);
        line_props(plotH, sig)
    end
    pretty_plot(axisH)
    xlabel("time [sec]",FontSize=15, Interpreter="none")
    title("Time Analysis");
    ylim(timeLimits)

    %% Plot freqs:
    [L, R, f] = Algo.filteredFFTs(left, right);    
    axisH = subplot(1,3,2);
    cellItems = {{L,left},{R,right}};
    for iC = 1 : length(cellItems)
        X = cellItems{iC}{1};
        sig = cellItems{iC}{2};
        hold(axisH, "on");
        plotH = plot(f, abs(X));
        line_props(plotH, sig)
    end
    pretty_plot(axisH)
    xlabel("freq [Hz]", FontSize=15, Interpreter="none")
    ylabel("|fft|"    , FontSize=15, Interpreter="none")
    title("Frequency Analysis");
    ylim(freqLimits)

    %% Plot BSI:
    axisH = subplot(1,3,3);
    bsi = Algo.singleBSI(L,R);
    plotH = bar(0, bsi,1, "stacked", "blue");
    xlim([0,1/2]);
    ylim([-1,1]);
    title("BSI")
    %% Adjust sizes:
    subaxes = get_subaxes(figH);    
    widen_axes(subaxes, 1, 1.5);
    widen_axes(subaxes, 2, 1.5);
    widen_axes(subaxes, 3, 0.1);
    figH.Position(1) = figH.Position(1)-200;
    figH.Position(2) = figH.Position(2)-200;
    
    %% Capture video:
    figRecorder.capture(figH);
    %% Update graphics:
    close(figH);
    progBar.step()
end

%% Close:
progBar.close();
figRecorder.finish();
disp("End.")
%% End:

%% Subs:
function line_props(plotH, sig) 
    side  = string(sig.placement.side);
    label = string(sig.placement.name);
    plotH.DisplayName = label + " (" + side + ")";
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

        % Derive current limits:
        crntLims.time.max = max([left.recording; right.recording]);
        crntLims.time.min = min([left.recording; right.recording]);
        crntLims.freq.max = max([abs(L); abs(R)]);
        crntLims.freq.min = min([abs(L); abs(R)]);
        crntLims.bsi.max  = bsi;
        crntLims.bsi.min  = bsi;

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
function subaxes = get_subaxes(figH)
    subaxes = [];
    children = figH.Children;
    for iC = 1 : length(children)
        child = children(iC);
        if string(child.Type) == "axes"
            subaxes = [child, subaxes];
        end
    end % for iC
end
%% 
function widen_axes(subaxes, ind, factor)
    % calc width:
    prevWidth = subaxes(ind).Position(3);
    newWidth = prevWidth * factor;
    % Assign:
    subaxes(ind).Position(3) = newWidth;
    % adjust next axes:
    dW = newWidth-prevWidth;
    for nextInd = ind+1 : length(subaxes)
        drawnow;
        subaxes(nextInd).Position(1) = subaxes(nextInd).Position(1) + dW;
    end
    
end