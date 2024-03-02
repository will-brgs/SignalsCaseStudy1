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

disp("leela says hi");
disp('Will says hi')
%% Task 1: Design Variable Amplification for 5 Band Frequencies
% Here we will construct the set of linear systems
%% Task 1: Initialize R and C values for desired cutoff frequencies
C = 10e-6; % Consider resistance will be constant at 10uF, R will change to alter cutoff freq
% Create a vector of R values for Lo and Hi respectively
R_Lo = zeros(5,1);
R_Hi = zeros(5,1);

center_band = [60, 230, 910, 3e3, 14e3]; % 5 band frequency centerpoints
cutoffs = zeros(length(center_band),2); % 1st column lo, 2nd hi


k_cut = 0.1; % Threshod of what magnitude a center frequency will extend past itself

%% R value Allocaton
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
clear i, clear k_cut; clear cutoff_Hi, clear cutoff_Lo
%% Lsim Coefficient Allocation
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

clear C;

%% Bode Plot Test: Independent Bandpass Filters
% Generate Bode magnitude plots for all 5 bandpass filters
bode_size = 200; % How many differnt frequencies we want to test
bode_freq = logspace(1, 4.25, bode_size); %Generate different frequency vals, 10^4.25 max yeilds close to limit for human hearing
t = 0:1/Fs:1; %Sample timepoint vector

H = zeros(bode_size,1);

figure, hold on
sgtitle('Bode Plot Outputs for 5 Bandpass Filters')
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
semilogx(bode_freq, H_mag, 'linewidth', 1.5)
xlim([bode_freq(1),bode_freq(end)])
xlabel('Frequency (Hz)');
ylabel('Output (dB)');
xline(center_band(j), "-");
xline(cutoffs(j,1), "--"); % Lower Frequency Cutoff
xline(cutoffs(j,2), "--"); % Upper frequency Cutoff
title("Bode Plot Magnitude",num2str(center_band(j)) + " Hz Center")
end
clear cutoff_Hi, clear cutoff_Lo, clear x_filter, clear i, clear j, clear x
clear H, clear H_mag,
%% Bode Plot Test 2: Combined Equalizer
% t, bode_freq, are resued from revious bode plots

H = zeros(bode_size,1);
x_sum = zeros(length(t), 1);

for i = 1:length(bode_freq)
    freq_current = bode_freq(i);
    x = exp(1j* 2*pi * freq_current * t);
    
    for j = 1:5
        x_out = lsim(b_Lo(j,:),a_Lo(j, :), x, t);
        x_out = lsim(b_Hi(j,:),a_Hi(j,:), x_out, t);
        x_sum = x_sum + x_out;
    end
    H(i) = x_sum(end)/x(end);
end

%Calculate magnitude and angle values of complex gain
H_mag = 20 * log(abs(H));

figure, hold on
semilogx(bode_freq, H_mag, 'linewidth', 1.5)
xlim([bode_freq(1),bode_freq(end)])
title("High-pass Magnitude")
for i = 1:5
    for j = 1:2
        xline(cutoffs(i,j), "--"); % Create cutoff lines for each center
    end
        xline(center_band(1,i), "-"); % Create centerlines
end
hold off
clear x, clear t, clear x_out, clear x_sum, clear H, clear H_mag
clear i, clear j, clear bode_size, clear bode_freq
%% Import Chirp Function
% Sample code from Hw2 to generate chirp between given frequencies
dT = 1/Fs; % sampling period
t = 0:dT:3; % time vector
fmin = 1; fmax = 23e3; % 10000; % min and max frequencies for chirp
chirp_f = (fmax-fmin).*t/max(t)+fmin; % chirp instantaneous frequency
chirp_x = cos(2*pi*chirp_f/2.*t); % chirp signal

clear dT, clear t,clear fmin, clear fmax, 
%% 
output = chirp_x;
outputSum = zeros(length(chirp_f), 1);

for i = 5:5
    output_filter = lsim(b_Lo(i,:),a_Lo(i, :), output, chirp_f);
    output_filter = lsim(b_Hi(i,:),a_Hi(i,:), output_filter, chirp_f);
    outputSum = outputSum + output_filter;
end

figure, plot(outputSum);
xscale log;
% for i = 1
% output_filter = lsim(b_Lo(i,:),a_Lo(i,:), output,fchirp);
% output_filter = lsim(b_Hi(i,:),a_Hi(i,:), output_filter,fchirp);
% outputSum = outputSum + output_filter;
% end

%% Read All Audio Files

%Import Blue in Green with Siren
[sound_BGS] = audioread('Blue in Green with Siren.wav','native');
sound_BGS = sound_BGS(:,1);

%Import Giant Steps Bass Cut
[sound_GSBC] = audioread('Giant Steps Bass Cut.wav','native');
sound_GSBC = sound_GSBC(:,1);

%Import piano_noisy
[sound_PN] = audioread('piano_noisy.wav','native');
sound_PN = sound_PN(:,1);

%Import roosevelt_noisy
[sound_RN] = audioread('roosevelt_noisy.wav','native');
sound_RN = sound_RN(:,1);

%Import Space Staion - Treble Cut
[sound_SS] = audioread('Space Station - Treble Cut.wav','native');
sound_SS = sound_SS(:,1);

%Import violin_w_siren
[sound_VS] = audioread('violin_w_siren.wav','native');
sound_VS = sound_VS(:,1);