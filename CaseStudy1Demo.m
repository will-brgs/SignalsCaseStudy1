%% Case Study 1 - Demo

%% load noisy audio data, then play recording
[xv,xvfs] = audioread('violin_w_siren.wav'); 
sound(xv,xvfs)

%% show time-domain waveform
t = [0:length(xv)-1]*1/xvfs; 
figure, plot(t,xv)
xlabel('t, seconds')
ylabel('amplitude')

%% show transform 
f = [0:length(xv)-1]*xvfs/length(xv);
XV = fft(xv); 
figure, plot(f,abs(XV));
xlabel('f, Hz')
ylabel('|X(f)|')

%% show transform for first second
xvSnip = xv(1:xvfs); 
f = [0:length(xvSnip)-1]*xvfs/length(xvSnip);
XVSNIP = fft(xvSnip); 
figure, plot(f,abs(XVSNIP));
xlabel('f, Hz')
ylabel('|X(f)|')
pause
set(gca,'YScale','log')

%% show spectrogram
% figure, spectrogram(xv,256,200,256,xvfs)
figure, spectrogram(xv,1024,200,1024,xvfs)

%% second audio sample 
[xv,xvfs] = audioread('roosevelt_noisy.wav');
sound(xv,xvfs)
figure, spectrogram(xv,1024,200,1024,xvfs)

%% last audio sample 
[xv,xvfs] = audioread('piano_noisy.wav');
sound(xv,xvfs)
figure, spectrogram(xv,1024,200,1024,xvfs)


