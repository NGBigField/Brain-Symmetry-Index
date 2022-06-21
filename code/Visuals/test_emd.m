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

%% Bivariate Empirical Mode Decomposition
sigs = [left.recording, right.recording];
[IMFs, res ,info] = emd(left.recording);

%% Plot:
figH = figure;
sgtitle("EMD of left-brain EEG")
for i = 1 : 7
    imf = IMFs(:,i);
    label = "imf-"+string(i);
    subplot(4,2,i);
    p = plot(time, imf);
    p.DisplayName = label;    
    ylabel(label, FontSize=14)
end
subplot(4,2,8);
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