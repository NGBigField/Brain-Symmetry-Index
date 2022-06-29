% Before running this script, have one tagged-Data loaded to your workspace.
% These can be found in folder: Data\tagged\
%
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

%% Empirical Mode Decomposition
[IMFs, res ,info] = emd(left.recording);

%% Plot:
figH = figure;
m=5; n=2;
sgtitle("EMD of left-brain EEG")
subplot(m,n,1)
plot(time, left.recording);
ylabel("Input", FontSize=14)
for i = 1 : 7
    imf = IMFs(:,i);
    label = "imf-"+string(i);
    subplot(m,n,2+i);
    p = plot(time, imf);
    p.DisplayName = label;    
    ylabel(label, FontSize=14)
end
subplot(m,n,10);
p = plot(time, res);
ylabel("Res", FontSize=14)

%% End

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