%% Signals and Systems Case Study 1: Audio Equalizer
%% Introduction
% * Authors:                  Will Burgess, Leela Srinivas, Mack La Rosa
% * Class:                    ESE 351
% * Date:                     Created 2/20/2024, Last Edited 
% * With contributions from:  
%% Initialization
close all
Fs = 44.1e3; %44.1 kHz Audio Sampling Frequency

disp("leela says hi");

%% Task 1: Design Variable Amplification for 5 Band Frequencies
% Here we will construct the set of linear systems
%% Task 1: Initialize R and C values for desired cutoff frequencies
C = 10e-6; % Consider resistance will be constant at 10uF, R will change to alter cutoff freq
% Create a vector of R values for Lo and Hi respectively
R_Hi = zeros(5,1);
R_Lo = zeros(5,1);
% Band 1: 60Hz
cutoff_Lo = 1;
R_Lo(1) = 1/(2*pi*C*cutoff_Lo);

cutoff_Hi = 119;
R_Hi(1) = 1/(2*pi*C*cutoff_Hi);

%Band 2: 230Hz
cutoff_Lo = 119;
R_Lo(2) = 1/(2*pi*C*cutoff_Lo);

cutoff_Hi = 341;
R_Hi(2) = 1/(2*pi*C*cutoff_Hi);

%Band 3: 910Hz
cutoff_Lo = 341;
R_Lo(3) = 1/(2*pi*C*cutoff_Lo);

cutoff_Hi = 1479;
R_Hi(3) = 1/(2*pi*C*cutoff_Hi);

%Band 4: 3kHz
cutoff_Lo = 1479;
R_Lo(4) = 1/(2*pi*C*cutoff_Lo);

cutoff_Hi = 4521;
R_Hi(4) = 1/(2*pi*C*cutoff_Hi);

%band 5: 14kHz
cutoff_Lo = 4521;
R_Lo(5) = 1/(2*pi*C*cutoff_Lo);

cutoff_Hi = 23479;
R_Hi(5) = 1/(2*pi*C*cutoff_Hi);

%% Task 1: Initialize a and b Coefficients for lsim of bands

%Lowpass Filter Coefficients
a_Lo = zeros(5, 2);
a_Lo(:,1) = 1;
a_Lo(:,2) = 1./(C.*R_Lo);
b_Lo = 1./(C.*R_Lo);

%HighPass Filter Coefficients
a_Hi = zeros(5, 1);
a_Hi(:,1) = 1;
a_Hi(:,2) = 1./(C.*R_Lo);
b_Hi = zeros(5, 1);
b_Hi(:,1) = 1;

%% Import Chirp

fs = 44.1e3; % sampling frequency
dT = 1/fs; % sampling period
t = 0:dT:3; % time vector
fmin = 1; fmax = 23e3; % 10000; % min and max frequencies for chirp
fchirp = (fmax-fmin).*t/max(t)+fmin; % chirp instantaneous frequency
xchirp = cos(2*pi*fchirp/2.*t); % chirp signal
% xchirp = chirp(t,fmin,max(t),fmax,'logarithmic'); % use of Matlab chirp with logarithmic frequency variation
% sound(xchirp,fs)
ychirpHPF = filter(.5*[1 -1],1,xchirp);
ychirpMA = filter(ones(2,1),1,xchirp);
% visualization - log-frequency linear-amplitude
figure, subplot(3,1,1), plot(fchirp,xchirp); title('Chirp input')
subplot(3,1,2), plot(fchirp,(ychirpMA)); title('DT Moving Average')
subplot(3,1,3), plot(fchirp,(ychirpHPF)); title('DT Highpass Filter'), xlabel('Frequency (Hz)')
% to visualize with log or linear scale for frequency or amplitude, use this code accordingly
for i = 1:3, subplot(3,1,i), set(gca,'XScale','log'), set(gca,'YScale','linear'), axis tight, end
%% Task 1: 
output = xchirp;
for i = 1:5
output = lsim(b_Lo(i,:),a_Lo(i,:), output,fchirp);
output = lsim(b_Hi(i,:),a_Hi(i,:), output,fchirp);
end
figure, plot(output)
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