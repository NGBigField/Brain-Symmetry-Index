%% clear
clearvars('-except','Data');
close all; clc;

%% Prepare:
figRecorder = Classes.FigureRecorder(); % Prepare Recorder
% Adjust plot limits:
[timeLimits, freqLimits] = derive_plot_limits(Data.EEGs);
timeLimits = timeLimits/2;
freqLimits = freqLimits/4;
freqLimits(1) = 0;

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
    axisH = subplot(1,2,1);
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
    axisH = subplot(1,2,2);
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
function [timeLimits, freqLimits] = derive_plot_limits(EEGs)
    %% Init limits:
    timeLimits = [inf, -inf];
    freqLimits = [inf, -inf];

    %% Iterate:
    nData = length(EEGs);
    for i = 1 : nData
        % get data:
        eegData = EEGs(i);
        left  = eegData.left;
        right = eegData.right;
        [L, R, ~] = Algo.filteredFFTs(left, right);      
        % Derive current limits:
        tMax = max([left.recording; right.recording]);
        tMin = min([left.recording; right.recording]);
        fMax = max([abs(L); abs(R)]);
        fMin = min([abs(L); abs(R)]);

        % Derive global limits:
        timeLimits(1) = min( timeLimits(1), tMin);
        timeLimits(2) = max( timeLimits(2), tMax);
        freqLimits(1) = min( freqLimits(1), fMin);
        freqLimits(2) = max( freqLimits(2), fMax);
    end
end