%% clear
clearvars('-except','Data');
addpath(genpath(pwd));
close all; clc;

%% Constants:
side_labels = ["left", "right"];
pre_filter = true;

%% parse data:
eeg_data = Data.EEGs(101);
left  = eeg_data.left;
right = eeg_data.right;
time  = eeg_data.time;   
fs    = eeg_data.left.props.frequency;

%% Filter signals:
if pre_filter
    l_sig = DSP.band_pass_filter(left.recording , fs);
    r_sig = DSP.band_pass_filter(right.recording, fs);
else
    l_sig = left.recording;
    r_sig = right.recording;
end

%% Plot filter results:
if pre_filter
    plot_filter_results(left.recording, l_sig, fs)
end

%% Bivariate Empirical Mode Decomposition
num_components = 8;
X = [ l_sig, r_sig] ;
[imf_tensor, residue] = Algo.bemd(X, num_components=num_components);
num_components = size(imf_tensor, 3);

%% Plot IMFs:
figH = figure;
sgtitle("BEMD EEG"+newline+"condition="+string(eeg_data.condition)+newline+filter_string(pre_filter))
for i = 0 : (num_components+1)
    % IMFS:
    for j = 1 : 2
        side_label = side_labels(j);       

        if i == 0 % Base
            sig = X(:,j);
            label = "Input";
        elseif i == num_components + 1
            sig = residue(:,j);
            label = "Res";
        else % IMFs
            sig = imf_tensor(:,j,i);
            label = "imf-"+string(i);
        end        

        subplot(num_components+2, 2, i*2 + j );
        plot(time, sig);

        % labels:
        if j == 1
            ylabel(label, FontSize=14)
        end
        if i == 0
            title(side_label, FontSize=18)
        end
    end % for j    

     
end % for i
% Shift size:
Visuals.standart_bemd_fig_size(figH)
% Connect x axes:
Visuals.link_axes(figH, "x")

%% Plot HHTs:
figH = figure;
sgtitle("BEMD EEG Hilbert-Huang Transform "+newline+"condition="+string(eeg_data.condition)+newline+filter_string(pre_filter))
color_oredr = colororder;

for iSide = 1 : 2
    axisH = subplot(1,2,iSide);
    title( side_labels(iSide) + " IMFs", FontSize=14)

    for iImf = 1 : num_components        
        color = color_oredr(iImf, :);           
        imf = imf_tensor(:, iSide, iImf);
        name = "imf-"+string(iImf);
        plot_hht(axisH, fs, imf, color, name)
        draw_pause()

    end % imfF

end % iSide
Visuals.link_axes(figH, "xy");
legend(Location="northwest")

%% Save figs:
%{
    Visuals.save_all_figs()
%}

%% Finish:
disp("Done.")
Sounds.finish()

%% End

%% Subs:

function plot_filter_results(raw, filtered, fs)
    % time vector:
    Ts = 1/fs;
    N = length(raw);
    t_vec = 0 : Ts : (N-1)*Ts;

    % time plot:
    figure;
    subplot(1,2,1)
    hold on; plot(t_vec, raw,      DisplayName="Raw"     )
    hold on; plot(t_vec, filtered, DisplayName="Filtered")
    xlabel("time [sec]", FontSize=16)
    ylabel("$$ x(t) $$", Interpreter='latex', FontSize=16)
    title("Time Domain", FontSize=16)
    legH = legend();
    legH.FontSize = 14;

    % freq plot:
    subplot(1,2,2)
    [Raw, f_raw] = Algo.fft(raw, fs);
    [Filt, f_filt] = Algo.fft(filtered, fs);
    hold on; plot(f_raw,  abs(Raw),  DisplayName="Raw"     )
    hold on; plot(f_filt, abs(Filt), DisplayName="Filtered")
    title("Frequency Domain", FontSize=16)
    xlabel(" f [Hz] ", FontSize=16)
    ylabel("$$ | X^{F}(f) | $$", Interpreter='latex', FontSize=16)
    legH = legend();
    legH.FontSize = 14;
end

%%

function plot_hht(axisH, fs, imf, color, name)
    [t, f, h] = Algo.hht(imf, fs);
    color_vals = Visuals.color_vector(h, color=color);

    % Plot:
    hold(axisH, "on")
    scatter(axisH, t, f, 3, color_vals, "filled", DisplayName=name)
    xlabel("time [sec]", FontSize=16)
    ylabel("freq [Hz]" , FontSize=16)
end

%%

function draw_pause()
    drawnow();
    pause(0.001);
end

%%

function str = filter_string(filter_on)
    if filter_on
        str = "with band-pass filter";
    else
        str = "without filter";
    end
end