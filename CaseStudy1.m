%% Signals and Systems Case Study 1: Audio Equalizer
%% Introduction
% * Authors:                  Will Burgess, Leela Srinivas, Mack La Rosa
% * Class:                    ESE 351
% * Date:                     Created 2/20/2024, Last Edited 
% * With contributions from:  
%% Initialization
close all
clear
Fs = 44.1e3; %44.1 kHz Audio Sampling Frequency
font_size = 14;
%% Design Task 1: Design Variable Amplification for 5 Band Frequencies
% The nessecary framework for the equalizer is constructed here
%% Design Task 1: Initialize R and C values for desired cutoff frequencies
C = 10e-6; % Consider resistance will be constant at 10uF, R will change to alter cutoff freq
% Create a vector of R values for Lo and Hi respectively
R_Lo = zeros(5,1);
R_Hi = zeros(5,1);

center_band = [60, 230, 910, 3e3, 14e3]; % 5 band frequency centerpoints
cutoffs = zeros(length(center_band),2); % 1st column lo, 2nd hi

k_cut = 0.25; % Threshod of what magnitude a center frequency will extend past itself

%% Design Task 1: R value Allocaton
for i = 1:5 % Calculate and allocate A and B values for Hi and Lo Coefficients

% Calc Hipass R
cutoff_Hi = center_band(i) - (k_cut * center_band(i));
cutoffs(i,1) = cutoff_Hi;
R_Hi(i) = 1/(2 * pi * C * cutoff_Hi);

% Calc Lopass R
cutoff_Lo = center_band(i) + (k_cut * center_band(i));
cutoffs(i,2) = cutoff_Lo;
R_Lo(i) = 1/(2 * pi * C * cutoff_Lo);

end
%% Design Task 1: Lsim Coefficient Allocation
%Lowpass Filter Coefficients
a_Lo = zeros(5, 2);
a_Lo(:,1) = 1;
a_Lo(:,2) = 1./(C.*R_Lo);

b_Lo = 1./(C.*R_Lo);

%HighPass Filter Coefficients
a_Hi = zeros(5, 1);
a_Hi(:,1) = 1;
a_Hi(:,2) = 1./(C.*R_Hi);

b_Hi = zeros(5, 2);
b_Hi(:,1) = 1;

%% Design Task 2: Test lsim variables by frequency response analysis
% Various bode plots are constructed to quantize the frequency response of
% the equalizer system
%% Design Task 2a: Bode Plot Test: Independent Bandpass Filters
% Generate Bode magnitude plots for all 5 bandpass filters
bode_size = 100; % How many differnt frequencies we want to test
bode_freq = logspace(1, 4.25, bode_size); %Generate different frequency vals, 10^4.25 max yeilds close to limit for human hearing
t = 0:1/Fs:0.25; %Sample timepoint vector

H = zeros(bode_size,1);

% Create bode plot for magnitude response
figure, hold on
sgtitle('Bode Plot Magnitude Response for 5 Bandpass Filters')
for j = 1:length(center_band)
for i = 1:length(bode_freq) % Generate bode plot
    freq_current = 2*pi * bode_freq(i); %Convert logspace freq to angular
    x = exp(1j * freq_current * t);

    % Pass through Both hi and lowpass filters to create bandpass
    x_filter = lsim(b_Lo(j,:),a_Lo(j,:), x, t);
    x_filter = lsim(b_Hi(j,:),a_Hi(j,:), x_filter, t);

    H(i) = x_filter(end)/x(end);
end

%Calculate magnitude and angle values of complex gain
H_mag = 20 * log(abs(H)); % Convert output H to dB magnitude

subplot(3,2,j)
semilogx(bode_freq, H_mag, 'linewidth', 2.25)
xlim([bode_freq(1),bode_freq(end)])
xlabel('Frequency (Hz)');
ylabel('Output (dB)');
xline(center_band(j), "-");
xline(cutoffs(j,1), "--"); % Lower Frequency Cutoff
xline(cutoffs(j,2), "--"); % Upper frequency Cutoff
title("Magnitude Response",num2str(center_band(j)) + " Hz Center Frequency")
end

% Repeat plotting but for phase response,
figure, hold on
sgtitle('Bode Plot Phase Response for 5 Bandpass Filters')
for j = 1:length(center_band)
for i = 1:length(bode_freq) % Generate bode plot
    freq_current = 2*pi * bode_freq(i); %Convert logspace freq to angular
    x = exp(1j * freq_current * t);

    % Pass through Both hi and lowpass filters to create bandpass
    x_filter = lsim(b_Lo(j,:),a_Lo(j,:), x, t);
    x_filter = lsim(b_Hi(j,:),a_Hi(j,:), x_filter, t);

    H(i) = x_filter(end)/x(end);
end

%Calculate angle values of complex gain
H_phase =  angle(H); % Convert output H to dB magnitude

subplot(3,2,j)
semilogx(bode_freq, H_phase, 'linewidth', 2.25, 'color','r')
xlim([bode_freq(1),bode_freq(end)])
xlabel('Frequency (Hz)');
ylabel('Phase (Radians)');
xline(center_band(j), "-", 'LineWidth', 1.5);
xline(cutoffs(j,1), "--", 'LineWidth', 1.5); % Lower Frequency Cutoff
xline(cutoffs(j,2), "--", 'LineWidth', 1.5); % Upper frequency Cutoff
title("Bode Plot Phase Response",num2str(center_band(j)) + " Hz Center Frequency")
end
%% Design Task 2b: Bode Plot Test 2: Combined Equalizer
% t, bode_freq, are reused from previous bode plots
H = zeros(bode_size,1);
gains = [5 2 5 2 5] * 1/4;
t = 0:1/Fs:0.25;

for i = 1:length(bode_freq)
    freq_current = bode_freq(i);
    x = exp(1j* 2*pi * freq_current * t);
    x_sum = zeros(length(t), 1);
    
        for j = 1:5
            x_out = lsim(b_Lo(j,:),a_Lo(j, :), x, t);
            x_out = lsim(b_Hi(j,:),a_Hi(j,:), x_out, t);
            x_sum = x_sum + gains(j) * x_out;
        end
        
    H(i) = x_sum(end)/x(end);
end

%Calculate magnitude and angle values of complex gain
H_mag = 20 * log(abs(H));
H_phase = (angle(H));

figure, hold on
subplot(2,1,1)
plot(bode_freq, H_mag, 'linewidth', 2.25,'color','b') % For some reason semilogx doesnt work here
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)', 'FontSize', font_size);
ylabel('Output (dB)', 'FontSize', font_size);
xlim([bode_freq(1),bode_freq(end)])
title("Merged Bandpass Equalizer Magnitude Output", 'FontSize', font_size); 
%legend(num2str(filter_n_times) + " times", num2str(filter_m_times) + " times");
for i = 1:5
    for j = 1:2
        xline(cutoffs(i,j), "--", 'LineWidth', 1.5); % Create cutoff lines for each center
    end
        xline(center_band(1,i), "-", 'LineWidth', 1.5); % Create centerlines
end

subplot(2,1,2)
plot(bode_freq, H_phase, 'linewidth', 2.25,'color','r') 
set(gca, 'XScale', 'log');
font_size = 14;
xlabel('Frequency (Hz)', 'FontSize', font_size);
ylabel('Phase (Radians)', 'FontSize', font_size);
xlim([bode_freq(1),bode_freq(end)])
title("Merged Bandpass Equalizer Phase Output", 'FontSize', font_size); 
%legend(num2str(filter_n_times) + " times", num2str(filter_m_times) + " times");
for i = 1:5
    for j = 1:2
        xline(cutoffs(i,j), "--", 'LineWidth', 1.5); % Create cutoff lines for each center
    end
        xline(center_band(1,i), "-", 'LineWidth', 1.5); % Create centerlines
end
hold off
%% Design Task 2b: Merged Bode Plot Test, ends pulled up(unsucessful)
% t, bode_freq, are resued from revious bode plots
H = zeros(bode_size,1);
%filter_m_times = 5;
gains = [1 1 1 1 1];
t = 0:1/Fs:0.25;

for i = 1:length(bode_freq)
    freq_current = bode_freq(i);
    x = exp(1j* 2*pi * freq_current * t);
    x_sum = zeros(length(t), 1);
    
    for j = 1:5
        if j == 1
        %lsim low
        x_band = lsim(b_Lo(j,:),a_Lo(j, :), x, t);

        elseif j == 5
        %Lsim hi
        x_band = lsim(b_Hi(j,:),a_Hi(j,:), x, t);

        else % ie i = 2:4
        %lsim low
        x_lo = lsim(b_Lo(j,:),a_Lo(j, :), x, t);
        %Lsim hi
        x_band = lsim(b_Hi(j,:),a_Hi(j,:), x_lo, t);
        end
    
            x_sum = x_sum + gains(j)*x_band;
    end      
    H(i) = x_sum(end)/x(end);
end

%Calculate magnitude and angle values of complex gain
H_mag = 20 * log(abs(H));
%H_mag = H_mag(1:end); %Remove initial outlier data point(outside of human hearing)

figure, hold on
subplot(2, 1, 1)
plot(bode_freq, H_mag, 'linewidth', 2.25,'color','b') 
set(gca, 'XScale', 'log');
font_size = 14;
xlabel('Frequency (Hz)', 'FontSize', font_size);
ylabel('Output (dB)', 'FontSize', font_size);
xlim([bode_freq(1),bode_freq(end)])
title("Ends Pulled Up Bandpass Equalizer Magnitude", 'FontSize', font_size); 
%legend(num2str(filter_n_times) + " times", num2str(filter_m_times) + " times");
for i = 1:5
    for j = 1:2
        xline(cutoffs(i,j), "--", 'LineWidth', 1.5); % Create cutoff lines for each center
    end
        xline(center_band(1,i), "-", 'LineWidth', 1.5); % Create centerlines
end

subplot(2, 1, 2)
plot(bode_freq, H_phase, 'linewidth', 2.25,'color','r') 
set(gca, 'XScale', 'log');
font_size = 14;
xlabel('Frequency (Hz)', 'FontSize', font_size);
ylabel('Phase (Radians)', 'FontSize', font_size);
xlim([bode_freq(1),bode_freq(end)])
title("Ends Pulled Up Bandpass Equalizer Phase", 'FontSize', font_size); 
%legend(num2str(filter_n_times) + " times", num2str(filter_m_times) + " times");
for i = 1:5
if i == 1
    xline(cutoffs(i,2), "--", 'LineWidth', 1.5); % Create cutoff only on high side
elseif i == 5
    xline(cutoffs(i,1), "--", 'LineWidth', 1.5); % Create cutoff only on low side
else % ie the three center lines
    for j = 1:2
    xline(cutoffs(i,j), "--", 'LineWidth', 1.5); % Create cutoff on both sides
    end
end
xline(center_band(1,i), "-", 'LineWidth', 1.5); % Create centerlines
end
hold off
%NOTE: From herein out a function will replace the implementation of the
%equalizer linear simulation, this function is denoted EqualizerFunc
%% Testing Task 1: Design 3 Audio Presets
% Here gain presets are chosen to form the merged bandpass equalizer
% frequency repsonse. This analysis is done via bode plots. Note that phase
% response is not calculated as it is irrelevent to the way the audio is
% desired to be filtered
%% Testing Task 1: Treble Boost
% t, bode_freq, are resued from revious bode plots
H = zeros(bode_size,1);
gains_treble = [1 1 0.75 3 4];
t = 0:1/Fs:0.25;

for i = 1:length(bode_freq)
    freq_current = bode_freq(i);
    x = exp(1j* 2*pi * freq_current * t);
    x_sum = zeros(length(t), 1);
    
    for j = 1:5
        x_out = lsim(b_Lo(j,:),a_Lo(j, :), x, t);
        x_out = lsim(b_Hi(j,:),a_Hi(j,:), x_out, t);
        x_sum = x_sum + gains_treble(j) * x_out;
    end
    H(i) = x_sum(end)/x(end);
end

%Calculate magnitude and angle values of complex gain
H_mag = 20 * log(abs(H));
H_phase = (angle(H));

figure, hold on
subplot(2, 1, 1)
plot(bode_freq, H_mag, 'linewidth', 2.25) % For some reason semilogx doesnt work here
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)', 'FontSize', font_size);
ylabel('Output (dB)', 'FontSize', font_size);
xlim([bode_freq(1),bode_freq(end)])
title("Treble Boost Equalizer Magnitude", 'FontSize', font_size); 
%legend(num2str(filter_n_times) + " times", num2str(filter_m_times) + " times");
for i = 1:5
    for j = 1:2
        xline(cutoffs(i,j), "--", 'LineWidth', 1.5); % Create cutoff lines for each center
    end
        xline(center_band(1,i), "-", 'LineWidth', 1.5); % Create centerlines
end


subplot(2, 1, 2)
plot(bode_freq, H_phase, 'linewidth', 2.25, 'Color', 'r') % For some reason semilogx doesnt work here
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)', 'FontSize', font_size);
ylabel('Output (dB)', 'FontSize', font_size);
xlim([bode_freq(1),bode_freq(end)])
title("Treble Boost Equalizer Phase", 'FontSize', font_size);
for i = 1:5
    for j = 1:2
        xline(cutoffs(i,j), "--", 'LineWidth', 1.5); % Create cutoff lines for each center
    end
        xline(center_band(1,i), "-", 'LineWidth', 1.5); % Create centerlines
end
hold off

%% Testing Task 1: Bass Boost
% t, bode_freq, are resued from revious bode plots
H = zeros(bode_size,1);
gains_bass = [4 3 0.75 1 1.5];
t = 0:1/Fs:0.25;

for i = 1:length(bode_freq)
    freq_current = bode_freq(i);
    x = exp(1j* 2*pi * freq_current * t);
    x_sum = zeros(length(t), 1);
    
    for j = 1:5
        x_out = lsim(b_Lo(j,:),a_Lo(j, :), x, t);
        x_out = lsim(b_Hi(j,:),a_Hi(j,:), x_out, t);
        x_sum = x_sum + gains_bass(j) * x_out;
    end
    H(i) = x_sum(end)/x(end);
end

%Calculate magnitude and angle values of complex gain
H_mag = 20 * log(abs(H));
H_phase = angle(H);

figure, hold on
subplot(2, 1, 1)
plot(bode_freq, H_mag, 'linewidth', 2.25) % For some reason semilogx doesnt work here
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)', 'FontSize', font_size);
ylabel('Output (dB)', 'FontSize', font_size);
xlim([bode_freq(1),bode_freq(end)])
title("Bass Boost Equalizer Magnitude", 'FontSize', font_size); 
%legend(num2str(filter_n_times) + " times", num2str(filter_m_times) + " times");
for i = 1:5
    for j = 1:2
        xline(cutoffs(i,j), "--", 'LineWidth', 1.5); % Create cutoff lines for each center
    end
        xline(center_band(1,i), "-", 'LineWidth', 1.5); % Create centerlines
end

subplot(2, 1, 2)
plot(bode_freq, H_phase, 'linewidth', 2.25, 'Color', 'r') % For some reason semilogx doesnt work here
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)', 'FontSize', font_size);
ylabel('Output (dB)', 'FontSize', font_size);
xlim([bode_freq(1),bode_freq(end)])
title("Bass Boost Equalizer Phase", 'FontSize', font_size);
for i = 1:5
    for j = 1:2
        xline(cutoffs(i,j), "--", 'LineWidth', 1.5); % Create cutoff lines for each center
    end
        xline(center_band(1,i), "-", 'LineWidth', 1.5); % Create centerlines
end
hold off

%% Testing Task 1: Unity
% t, bode_freq, are resued from revious bode plots
H = zeros(bode_size,1);
filter_n_times = 1;
%filter_m_times = 5;
gains_unity = [1 1 1 1 1];
t = 0:1/Fs:0.25;

for i = 1:length(bode_freq)
    freq_current = bode_freq(i);
    x = exp(1j* 2*pi * freq_current * t);
    x_sum = zeros(length(t), 1);
    
    for j = 1:5
        x_out = lsim(b_Lo(j,:),a_Lo(j, :), x, t);
        x_out = lsim(b_Hi(j,:),a_Hi(j,:), x_out, t);
        x_sum = x_sum + gains_unity(j) * x_out;
    end      
    H(i) = x_sum(end)/x(end);
end

%Calculate magnitude and angle values of complex gain
H_mag = 20 * log(abs(H));
H_mag = H_mag(2:end); %Remove initial outlier data point(outside of human hearing)
H_phase = angle(H);

figure, hold on
subplot(2,1,1)
plot(bode_freq(2:end), H_mag, 'linewidth', 2.25) % For some reason semilogx doesnt work here
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)', 'FontSize', font_size);
ylabel('Output (dB)', 'FontSize', font_size);
xlim([bode_freq(1),bode_freq(end)])
title("Unity Equalizer Magnitude", 'FontSize', font_size); 
%legend(num2str(filter_n_times) + " times", num2str(filter_m_times) + " times");
for i = 1:5
    for j = 1:2
        xline(cutoffs(i,j), "--", 'LineWidth', 1.5); % Create cutoff lines for each center
    end
        xline(center_band(1,i), "-", 'LineWidth', 1.5); % Create centerlines
end

subplot(2, 1, 2)
plot(bode_freq, H_phase, 'linewidth', 2.25, 'Color', 'r') % For some reason semilogx doesnt work here
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)', 'FontSize', font_size);
ylabel('Output (dB)', 'FontSize', font_size);
xlim([bode_freq(1),bode_freq(end)])
title("Unity Equalizer Phase", 'FontSize', font_size);
for i = 1:5
    for j = 1:2
        xline(cutoffs(i,j), "--", 'LineWidth', 1.5); % Create cutoff lines for each center
    end
        xline(center_band(1,i), "-", 'LineWidth', 1.5); % Create centerlines
end
hold off

%% Testing Task 2: Pass preset fitlers through Space Station and Giant Steps
% Here the same audio presets are passed through two sample clips. The code
% segnemnts are deisnged to be operated independently of one another.
%% Testing Task 2: Import Audio Files
%Import Space Station
[sound_SS,Fs_SS] = audioread('Space Station - Treble Cut.wav');
sound_SS = sound_SS(:,1);

%Import giant steps
[sound_GSBC,Fs_GSBC] = audioread('Giant Steps Bass Cut.wav');
sound_GSBC = sound_GSBC(:,1);
%% Testing Task 2: Space Station : Unity

out_SS = equalizerFunc(sound_SS, Fs_SS, gains_unity, center_band, k_cut);

sound(out_SS,Fs_SS);

figure, hold on
sgtitle('Spectrograms for Space Station', 'fontsize',font_size)

subplot(3,1,1)
[s,f,t] = spectrogram(out_SS,hamming(256),round(256/2),1024,Fs_SS);
imagesc(t, f, 20*log10(abs(s))); % Convert to dB scale for better visualization
axis xy; % Flip the y-axis to have low frequencies at the bottom
clim([-100,70]);
colorbar;
cbar = colorbar;
colormap default
ylabel(cbar,'Magnitude (dB)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Unity Spectrogram');
%% Testing Task 2: Space Station : Treble Boost

out_SS = equalizerFunc(sound_SS, Fs_SS, gains_treble, center_band, k_cut);

sound(out_SS,Fs_SS);
subplot(3,1,2)
s = spectrogram(out_SS,hamming(256),round(256/2),1024,Fs_SS);
imagesc(t, f, 20*log10(abs(s))); % Convert to dB scale for better visualization
axis xy; % Flip the y-axis to have low frequencies at the bottom
colormap default
colorbar;
cbar = colorbar;

clim([-100,70]);
ylabel(cbar,'Magnitude (dB)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Treble Boost Spectrogram');
%% Testing Task 2: Space Station : Bass Boost

out_SS = equalizerFunc(sound_SS, Fs_SS, gains_bass, center_band, k_cut);

sound(out_SS,Fs_SS);
subplot(3,1,3)
s = spectrogram(out_SS,hamming(256),round(256/2),1024,Fs_SS);
imagesc(t, f, 20*log10(abs(s))); % Convert to dB scale for better visualization
axis xy; % Flip the y-axis to have low frequencies at the bottom
colormap default
colorbar;
cbar = colorbar;
clim([-100,70]);
ylabel(cbar,'Magnitude (dB)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Bass Boost Spectrogram');
hold off
%% Testing Task 2: Giant Steps : Unity

out_GSBC = equalizerFunc(sound_GSBC, Fs_GSBC, gains_unity, center_band, k_cut);

sound(sound_GSBC,Fs_GSBC);
figure, hold on
sgtitle('Spectrograms for Giant Steps', 'fontsize',font_size)

subplot(3,1,1)
[s,f,t] = spectrogram(out_GSBC,hamming(256),round(256/2),1024,Fs_GSBC);
imagesc(t, f, 20*log10(abs(s))); % Convert to dB scale for better visualization
axis xy; % Flip the y-axis to have low frequencies at the bottom
colormap default
colorbar;
cbar = colorbar;
clim([-100,70]);
ylabel(cbar, 'Magnitude (dB)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Unity Spectrogram');
%% Testing Task 2: Giant Steps : Treble Boost

out_GSBC = equalizerFunc(sound_GSBC, Fs_GSBC, gains_treble, center_band, k_cut);

sound(sound_GSBC,Fs_GSBC);

subplot(3,1,2)
s = spectrogram(out_GSBC,hamming(256),round(256/2),1024,Fs_GSBC);
imagesc(t, f, 20*log10(abs(s))); % Convert to dB scale for better visualization
axis xy; % Flip the y-axis to have low frequencies at the bottom
colormap default
colorbar;
cbar = colorbar;
clim([-100,70]);
ylabel(cbar, 'Magnitude (dB)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Treble Boost Spectrogram');
%% Testing Task 2: Giant Steps : Bass Boost

out_GSBC = equalizerFunc(sound_GSBC, Fs_GSBC, gains_bass, center_band, k_cut);

sound(sound_GSBC,Fs_GSBC);

subplot(3,1,3)
s = spectrogram(out_GSBC,hamming(256),round(256/2),1024,Fs_SS);
imagesc(t, f, 20*log10(abs(s))); % Convert to dB scale for better visualization
axis xy; % Flip the y-axis to have low frequencies at the bottom
colormap default
colorbar;
cbar = colorbar;
clim([-100,70]);
ylabel(cbar, 'Magnitude (dB)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Bass Boost Spectrogram');
hold off
%% Testing Task 3: Filter out Background Noise
% A sample audio file of 'blue in green' has been provided, but with
% signficant background noise. The goal is to remove the background siren
% by creating a custom gain setting

%Import Blue in Green with Siren
[sound_BGS, Fs_BGS] = audioread('Blue in Green with Siren.wav');
sound_BGS = sound_BGS(:,1);

%% Testing Task 3: Create Spectrogram to Observe Problematic Frequencies
sound(sound_BGS,Fs_BGS);

figure, hold on
sgtitle('Spectrograms for Giant Steps w/ Siren', 'fontsize',font_size)

subplot(2,1,1)
[s,f,t] = spectrogram(sound_BGS,hamming(256),round(256/2),1024,Fs_BGS);
imagesc(t, f, 20*log10(abs(s))); % Convert to dB scale for better visualization
axis xy; % Flip the y-axis to have low frequencies at the bottom
colormap default
colorbar;
clim([-100,70]);
ylabel(cbar, 'Magnitude (dB)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Orignial Audio File Spectrogram');
%% Testing Task 3: Use Specialized Gain Filter to Remove Siren Noise
gains_BGS = [1 ,1e-2 ,1e-6 ,1e-2, 1];
center_band = [60, 730, 1546, 3e3, 14e3]; %adjust center band to attempt to cut out 1546Hz frequency
k_cut = 0.2;
out_BGS = equalizerFunc(sound_BGS, Fs_BGS, gains_BGS, center_band, k_cut);

sound(out_BGS,Fs_BGS);
%Generate spectrogram
subplot(2,1,2)
s = spectrogram(out_BGS,hamming(256),round(256/2),1024,Fs_BGS);
imagesc(t, f, 20*log10(abs(s))); % Convert to dB scale for better visualization
axis xy; % Flip the y-axis to have low frequencies at the bottom
colormap default
colorbar;
cbar = colorbar;
clim([-100,70]);
ylabel(cbar, 'Magnitude (dB)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Custom Filter Output Spectrogram');
hold off
%% Testing Task 4: Modify a Custom Sound File
% For this section a simple soundclip of somebody saying the word bruh was
% used. the goal was to recreate a trend of internet humor and increase the
% bass to massive levels. The goal was to break the audio as much as
% possible in the name of science.
%% Testing Task 4: Play original sound
[sound_Custom, Fs_Custom] = audioread('Bruh.mp4');
sound_Custom = sound_Custom(:,1);

sound(sound_Custom,Fs_Custom);
[s_og,f,t] = spectrogram(sound_Custom,hamming(256),round(256/2),1024,Fs_Custom);


%% Testing Task 4: Modify and play custom Output
gains_Custom = [1e4 ,0.001 ,0.001 ,2, 1];
center_band = [60, 230, 910, 3e3, 14e3];
k_cut = 0.2;
out_Custom = equalizerFunc(sound_Custom, Fs_BGS, gains_Custom, center_band, k_cut);
cbar = colorbar;
sound(out_Custom,Fs_Custom);
s_new = spectrogram(out_Custom,hamming(256),round(256/2),1024,Fs_Custom);

figure, hold on
sgtitle('Spectrograms for Custom Audio Modification', 'fontsize',font_size)
subplot(2,1,1)
imagesc(t, f, 20*log10(abs(s_og))); % Convert to dB scale for better visualization
axis xy; % Flip the y-axis to have low frequencies at the bottom
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Original Audio Clip Spectrogram');
colormap default
colorbar;
ylabel(cbar, 'Magnitude (dB)');
clim([-100,70]);

subplot(2,1,2)
imagesc(t, f, 20*log10(abs(s_new))); % Convert to dB scale for better visualization
axis xy; % Flip the y-axis to have low frequencies at the bottom
colormap default
colorbar;
cbar = colorbar;
ylabel(cbar, 'Magnitude (dB)');
clim([-100,70]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Modified Audio Clip Spectrogram');