% ===========================
% 1. Analysis - HFC, Spectral Flux, and Spectral Centroid(For Timbre) 
%
% 2. Alignment - Want to create an empty vector=to size of spec.
% Then HFC y(clean sample) and place the original y by its HFC location 
% after findpeaks into the empty vector at HFC locations of spec
%
% 3. Evaluation - HFC seems to be more exact for percussive instruments
% than what was proposed as the same but a more general approach,
% Spectral Flux [Bello, 2005].
% 
% 4. Comparison with the Existing System - Current data producing
% mistriggers at 1.92256(s) and 3.464069(s) in Steven Slate's Trigger 
% Software. The HFC or Spectral Flux algorithms I implemented with my peak 
% finding settings does not produce these errors. There were
% supposed to be 13 triggers - Slate's Trigger produced 15. On the same
% sample data, mine produced 13.
%
% Test on Snare Drum
% Nested for loop to place sample at locations
% Start Reding papers about drum transcription
% Upload to github
% ===========================

% ===========================
% ---- Load Sound File ----
% ===========================

wav = 'DANNY_THE TRIBAL_ACOUSTIC.wav';
mistrigger = 'DANNY_THE TRIBAL_MISTRIG.wav';
sample = 'DANNY_THE TRIBAL_CLEAN_SAMPLE.wav';
[x, fs] = audioread(wav);
 y = audioread(sample);
 z = audioread(mistrigger);

% ===========================
% ---- Initialize ----
% ===========================

T = .020;
N = T*fs;
nfft = 1024;
overlap = N/2;

% ===========================
% ---- Spectrogram ----
% ===========================

subplot(411)
specgram(x(:,1),nfft,fs,hamming(N),overlap);
spec=specgram(x(:,1),nfft,fs,hamming(N),overlap);
size(spec)
spectrogram(x(:,1),hamming(N),overlap,nfft,fs,'yaxis')
title(sprintf('Spectrogram %wav',overlap,size(spec,2)));

% ===========================
% ---- HFC ----
% ===========================

subplot(412)
HFC=zeros(1,size(spec,2));
for k = 2:size(spec,2)
    for j = 2 : size(spec,1)
        HFC(k) = HFC(k)+(abs(spec(j,k))).^2 * j;
    end
end

% Find peaks
Threshold = 20000;
Distance = .1;

ax=linspace(0,length(x)/fs,length(HFC));
plot(ax,HFC);
colorbar
xlabel('HFC')
axis([0,max(ax),0,max(HFC)])

[pks,locs] = findpeaks(HFC,ax,'MinPeakDistance',Distance,'MinPeakHeight',Threshold);
text(locs+.02,pks,num2str((1:numel(pks))'))

% ===========================
% ---- Spectral Flux ----
% ===========================

subplot(413)
spectral_flux=zeros(1,size(spec,2));
for k = 2:size(spec,2)
  currframe=spec(:,k);
  prevframe=spec(:,k-1);
  l2norm=sqrt(sum((abs(currframe)-abs(prevframe)).^2));
  spectral_flux(k)=l2norm;
end
ax=linspace(0,length(x)/fs,length(spectral_flux));
plot(ax,spectral_flux);
xlabel('Spectral Flux')
axis([0,max(ax),0,max(spectral_flux)])

% Find peaks
[pks,locs] = findpeaks(spectral_flux,ax,'MinPeakDistance',.1,'MinPeakHeight',75);
text(locs+.02,pks,num2str((1:numel(pks))'))

% ===========================
% ---- Spectral Centroid ----
% ===========================

subplot(414)
spectral_centroid=zeros(1,size(spec,2));
for k = 1:size(spec,2)
  frame = spec(:,k);
  spectral_centroid(k) = sum ((1:size(frame)).*abs(frame'));
end
ax=linspace(0,length(x)/fs,length(spectral_centroid));
plot(ax,spectral_centroid);
xlabel('Spectral Centroid')
axis([0,max(ax),0,max(spectral_centroid)])

% ===========================
% ---- Originals ----
% ===========================

figure(2)
subplot(311)
ax=linspace(0,length(x)/fs,length(x));
plot(ax,x)
xlabel('Original Acoustic')

figure(2)
subplot(312)
ax=linspace(0,length(z)/fs,length(z));
plot(ax,z);
xlabel('Mistrigger')

figure(2)
subplot(313)
ax=linspace(0,length(y)/fs,length(y));
plot(ax,y)
xlabel('Clean Sample')