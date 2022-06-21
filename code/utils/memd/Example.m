close all;
clear;
%% 
fs=1000; % should be > 10 times of the maximum frequency of interest
t=0:1/fs:10-1/fs;

f1=5;
f2=20;

X=sin(2*pi*f1*t);
Y=sin(2*pi*f1*t) + cos(2*pi*f2*t);
% X synchronises with Y at f1 (5Hz)

%% Perform NA-MEMD decomposition
No_noise_channels=10;

noise_process=zeros(length(X),No_noise_channels);

Noise_power=mean([var(X) var(Y)]);

for i=1:No_noise_channels
    noise_process(:,i)=randn(length(X),1).*sqrt(Noise_power);
end

memd_input=[X' Y' noise_process];

IMF=memd(memd_input);

%% Plot IMFs
column=3;
row=size(IMF,2);

figure('name','IMFs');
subplot(row+1,column,1);plot(t,memd_input(:,1));
title('X');
ylabel('Magnitude');
axis tight;
 ylim([-1 1]);

subplot(row+1,column,2);plot(t,memd_input(:,2));
title('Y');
ylabel('Magnitude');
axis tight;
ylim([-2 2]);

subplot(row+1,column,3);plot(t,memd_input(:,3));
title('WGN');
ylabel('Magnitude');
axis tight;
ylim([-1 1]);

for i=1:column
    for j=1:row

        subplot(row+1,column,((j)*column)+i);

        if i==1
            plot(t,squeeze(IMF(1,j,:)));
        elseif i==2
            plot(t,squeeze(IMF(2,j,:)));   
        elseif i==3
            plot(t,squeeze(IMF(3,j,:)));
        end

        ylabel(num2str(j));

        if j==size(IMF,2)
            xlabel('Time (s)');
        end

        axis tight;
        ylim([-1 1]);

    end
end

%% Plot PSDs
W_PSD=hamming(length(X));
NFFT=20000;

for i=1:size(IMF,2)
    [Pxx(i,:),f]=periodogram(squeeze(IMF(1,i,:)),W_PSD,NFFT,fs);
    [Pyy(i,:),f]=periodogram(squeeze(IMF(2,i,:)),W_PSD,NFFT,fs);
    [Pwgn(i,:),f]=periodogram(squeeze(IMF(3,i,:)),W_PSD,NFFT,fs);
end

figure;
subplot(3,1,1);
plot(f,10*log10(Pxx));
ylabel('PSD (dB/Hz)');
try
legend('IMF 1','IMF 2','IMF 3','IMF 4','IMF 5','IMF 6','IMF 7','IMF 8','IMF 9','IMF 10','IMF 11','IMF 12','IMF 13','IMF 14');
catch
    legend('IMF 1','IMF 2','IMF 3','IMF 4','IMF 5','IMF 6','IMF 7','IMF 8','IMF 9','IMF 10','IMF 11','IMF 12','IMF 13');
    try
        legend('IMF 1','IMF 2','IMF 3','IMF 4','IMF 5','IMF 6','IMF 7','IMF 8','IMF 9','IMF 10','IMF 11','IMF 12');
    catch
        legend('IMF 1','IMF 2','IMF 3','IMF 4','IMF 5','IMF 6','IMF 7','IMF 8','IMF 9','IMF 10','IMF 11');
    end
end
title('X');

subplot(3,1,2);
plot(f,10*log10(Pyy));
ylabel('PSD (dB/Hz)');
try
legend('IMF 1','IMF 2','IMF 3','IMF 4','IMF 5','IMF 6','IMF 7','IMF 8','IMF 9','IMF 10','IMF 11','IMF 12','IMF 13','IMF 14');
catch
    legend('IMF 1','IMF 2','IMF 3','IMF 4','IMF 5','IMF 6','IMF 7','IMF 8','IMF 9','IMF 10','IMF 11','IMF 12','IMF 13');
    try
        legend('IMF 1','IMF 2','IMF 3','IMF 4','IMF 5','IMF 6','IMF 7','IMF 8','IMF 9','IMF 10','IMF 11','IMF 12');
    catch
        legend('IMF 1','IMF 2','IMF 3','IMF 4','IMF 5','IMF 6','IMF 7','IMF 8','IMF 9','IMF 10','IMF 11');
    end
end
title('Y');

subplot(3,1,3);
plot(f,10*log10(Pwgn));
ylabel('PSD (dB/Hz)');
xlabel('Frequency (Hz)');
try
legend('IMF 1','IMF 2','IMF 3','IMF 4','IMF 5','IMF 6','IMF 7','IMF 8','IMF 9','IMF 10','IMF 11','IMF 12','IMF 13','IMF 14');
catch
    legend('IMF 1','IMF 2','IMF 3','IMF 4','IMF 5','IMF 6','IMF 7','IMF 8','IMF 9','IMF 10','IMF 11','IMF 12','IMF 13');
    try
        legend('IMF 1','IMF 2','IMF 3','IMF 4','IMF 5','IMF 6','IMF 7','IMF 8','IMF 9','IMF 10','IMF 11','IMF 12');
    catch
        legend('IMF 1','IMF 2','IMF 3','IMF 4','IMF 5','IMF 6','IMF 7','IMF 8','IMF 9','IMF 10','IMF 11');
    end
end
title('WGN');

%% Choose IMF with maximum power in channel 1
for i=1:size(IMF,2)
    Px(i)=var(squeeze(IMF(1,i,:)));
    Py(i)=var(squeeze(IMF(1,i,:)));
    Pwgn_2(i)=var(squeeze(IMF(1,i,:)));
end

[~,I_x]=sort(Px,'descend');
% [~,I_y]=sort(Py,'descend');
chosen_IMF=I_x(1);
% chosen_IMF=[I_y(2) I_y(1)];
 %% Calculate PSI
W_PSI=fs; %ideally should be at > 5 times of the longest period (lowest frequency)
W_OV_PSI=0.5;
 
temp=IMF(:,chosen_IMF,:);
% temp=IMF(:,[chosen_IMF(2) chosen_IMF(1)],:);
% temp=sum(temp,2);

[PSI,T_vec,Pairs]=MEMD_PSYNC_T(fs,W_PSI,W_OV_PSI,temp);
 
PSI_noise_ave=squeeze(mean(PSI(2:end,:,:)));
PSI_noise_std=squeeze(std(PSI(2:end,:,:)));

pair_XY=1;

figure;
plot(T_vec,squeeze(PSI(pair_XY,:,:)));
hold on;
plot(T_vec,PSI_noise_ave,'r');
hold on;
plot(T_vec,PSI_noise_ave+2*PSI_noise_std,'r--');
hold on;
plot(T_vec,PSI_noise_ave-2*PSI_noise_std,'r--');
xlabel('Time (s)');
ylabel('PSI');
axis tight;
ylim([0 1.1]);
legend('XY','Baseline','Upper bound','Lower bound');

%% Calculate TF representation of PSI 
FR_Hz = 0.1;
E_Th=0.0001;

% Do not change these 2 values
FR=FR_Hz/fs;
NFB=round(0.5/FR);

temp2=IMF(1:2,chosen_IMF,:);
% temp2=IMF(1:2,[chosen_IMF(2) chosen_IMF(1)],:);
% temp2=sum(temp,2);

% Since we here want to estimate synchrony between channels 1&2, the data channels, the first input argument of the 'MEMD_PSYNC_TF' function
% includes only the required channels (channels 1&2) -> memd_input(:,1,2).
 % The orignial data of the required channels is needed for determining threshold values
 [PSI_TF,F_vec,T_vec_2,Pairs]=MEMD_PSYNC_TF(memd_input(:,1:2),fs,W_PSI,FR,E_Th,NFB,W_OV_PSI,temp2);
 
figure;
imagesc(T_vec_2,fliplr(F_vec),flipud(squeeze(PSI_TF(1,:,:))));
set(gca,'YDir','normal');
colorbar;
colormap jet;
ylim([0 30]);
caxis([0 1])
xlabel('Time (s)');
ylabel('Frequency (Hz)');


 
 
 