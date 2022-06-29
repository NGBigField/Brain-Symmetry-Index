%% clear
clearvars('-except','Data');
close all; clc;

%% parse data:
eeg_data = Data.EEGs(1);
left  = eeg_data.left;
right = eeg_data.right;
time  = eeg_data.time;   

%% Plot sigs:
figH = figure;
axisH = axes(figH);
for sig = [left, right]
    record = sig.recording;
    hold(axisH, "on");
    plotH = plot(time, record);
    line_props(plotH, sig) 
end
pretty_plot(axisH)
xlabel("time [sec]")

%% Plot freqs:
[L, R, f] = Algo.filteredFFTs(left, right);
figH = figure;
axisH = axes(figH);
cellItems = {{L,left},{R,right}};
for iC = 1 : length(cellItems)
    X = cellItems{iC}{1};
    sig = cellItems{iC}{2};
    hold(axisH, "on");
    plotH = plot(f, abs(X));
    line_props(plotH, sig) 
end
pretty_plot(axisH)
xlabel("freq [Hz]")
ylabel("|fft|")
title("FFT elements");

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
    legend(axisH);
    grid(axisH,"on");
    grid(axisH,"minor");
end