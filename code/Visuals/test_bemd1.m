%% clear
clearvars('-except','Data');
addpath(genpath(pwd));
close all; clc;

%% Constants:
side_labels = ["left", "right"];
pre_filter = false;

%% parse data:
eeg_data = Data.EEGs(6);
left  = eeg_data.left;
right = eeg_data.right;
time  = eeg_data.time;   
fs    = eeg_data.left.props.frequency;

%% Filter signals:
if pre_filter
    f_pass1 = 2;
    f_pass2 = 45;
    high_pass_filt = @(x) highpass(x, f_pass1, fs);
    low_pass_filt = @(x) lowpass(x, f_pass2, fs);
    band_pass_filt = @(x) high_pass_filt(low_pass_filt(x));
    

    l_sig = band_pass_filt(left.recording);
    r_sig = band_pass_filt(right.recording);
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
sgtitle("BEMD EEG")
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

iImf = 3;
figH = figure;
sgtitle("BEMD EEG Hilbert-Huang Transform "+newline+"condition="+string(eeg_data.condition))

for iSide = 1 : 2
    % Basic plot data:
    axisH = subplot(1,2,iSide);
    title(side_labels(iSide) + "  imf-"+string(iImf), FontSize=15)
    drawnow();
    pause(0.001);

    % hht analysis:
    imf = imf_tensor(:, iSide, iImf);
    [hs, freq_used, time_used] = hht(imf, fs, FrequencyResolution=0.01);
    hs_max = full(max(hs, [], 'all'));
    hs_min = full(min(hs, [], 'all'));

    % Find elements:
    [non_zero_vec_i,  non_zero_vec_j, non_zero_hhs] = find(hs);
    f_vals = freq_used(non_zero_vec_i);
    t_vals = time_used(non_zero_vec_j);       
    color_limits = [0, 0.5];
    color_vals = Visuals.color_vector(non_zero_hhs, colormap="jet", value_limits=color_limits);

    % Plot:    
    hold(axisH, "on")
    scatter(axisH, t_vals, f_vals, 5, color_vals, "filled")
    xlabel("time [sec]", FontSize=16)
    ylabel("freq [Hz]" , FontSize=16)

    % Draw:
    drawnow();
    pause(0.001);
end
Visuals.link_axes(figH, "xy");

%%

%% Finish:
disp("Done.")
Sounds.gong(2, 110, 2)
Sounds.gong(2, 220, 3)

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


    

