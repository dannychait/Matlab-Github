% ===========================
%
% Danny Chait
% Onset Analysis and Resynthesis
% 
% ===========================

% ===========================
% ---- Load Sound File ----
% ===========================

wav = 'DANNY_THE TRIBAL_SNARE_ACOUSTIC.wav';
sample = 'DANNY_THE TRIBAL_SNARE_SAMPLE.wav';

[x, fs] = audioread(wav);
 y = audioread(sample);

% ===========================
% ---- Initialize ----
% ===========================

T = .020;
N = floor(T*fs);
nfft = 1024;
overlap = N/2;

% ===========================
% ---- Spectrograms ----
% ===========================

subplot(411)
specgram(x(:,1),nfft,fs,hamming(N),overlap);
spec=specgram(x(:,1),nfft,fs,hamming(N),overlap);
specSample=specgram(y(:,1),nfft,fs,hamming(N),overlap);
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
Threshold = 50000;
Distance = 1000;

ax=linspace(0,length(x),length(HFC));
plot(ax,HFC);
colorbar
xlabel('HFC')
axis([0,max(ax),0,max(HFC)])

[pks,locs] = findpeaks(HFC,ax,'MinPeakDistance',Distance,'MinPeakHeight',Threshold);
text(locs+.02,pks,num2str((1:numel(pks))'))

% ===========================
% ----- Alignment ----
% ===========================

% HFC Sample
HFCSample=zeros(1,size(specSample,2));
for k = 2:size(specSample,2)
    for j = 2 : size(specSample,1)
        HFCSample(k) = HFCSample(k)+(abs(specSample(j,k))).^2 * j;
    end
end

ax=linspace(0,length(x),length(HFCSample));
[pks1,locs1] = findpeaks(HFCSample,ax,'MinPeakDistance',Distance,'MinPeakHeight',Threshold);

figure(1)
subplot(413)
baseTime = 0;

% beat times in seconds
beatTimes = baseTime + (1/fs)*(locs'-1);

%attack_time = round(locs);
bufferLen = size(x,1);
sampleLen = numel(y);

sampleTrack = zeros(bufferLen+sampleLen,1);

% Place sample at each onset
for k = 1:numel(locs)
    idx = locs(k);
    sampleTrack(idx : idx+sampleLen-1) = y;
end

sampleTail = zeros(size(sampleLen,1),1);

% Overlap and add tails from previous frame
sampleTrack(1:sampleLen) = sampleTrack(1:sampleLen) + sampleTail;

ax=linspace(0,length(x),length(sampleTrack));
plot(sampleTrack)
colorbar
xlabel('Trigger Placement')
axis([min(ax),max(ax),min(sampleTrack),max(sampleTrack)])

% ===========================
% ---- Spectral Flux ----
% ===========================

figure(1)
subplot(414)
spectral_flux=zeros(1,size(spec,2));
for k = 2:size(spec,2)
  currframe=spec(:,k);
  prevframe=spec(:,k-1);
  l2norm=sqrt(sum((abs(currframe)-abs(prevframe)).^2));
  spectral_flux(k)=l2norm;
end
ax=linspace(0,length(x)/fs,length(spectral_flux));
plot(ax,spectral_flux);
colorbar
xlabel('Spectral Flux')
axis([0,max(ax),0,max(spectral_flux)])

[pks,locs] = findpeaks(spectral_flux,ax,'MinPeakDistance',.1,'MinPeakHeight',50);
text(locs+.02,pks,num2str((1:numel(pks))'))

% ===========================
% ---- Spectral Centroid ----
% ===========================

figure(3)
subplot(411)
spectral_centroid=zeros(1,size(spec,2));
for k = 1:size(spec,2)
  frame = spec(:,k);
  spectral_centroid(k) = sum ((1:size(frame)).*abs(frame'));
end
ax=linspace(0,length(x)/fs,length(spectral_centroid));
plot(ax,spectral_centroid);
colorbar
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
ax=linspace(0,length(y)*fs,length(y));
plot(ax,y)
xlabel('Clean Sample')