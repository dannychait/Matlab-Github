% ===========================
%
% Danny Chait
% Onset Analysis and Resynthesis
% 
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
N = floor(T*fs);
nfft = N;
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
Threshold = 80000;
Distance = 1000;

ax = 0:N-overlap:(N-overlap) * (length(HFC) - 1);
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

ax = 0:N-overlap:(N-overlap) * (length(HFCSample) - 1);
[pks1,locs1] = findpeaks(HFCSample,ax,'MinPeakDistance',Distance,'MinPeakHeight',Threshold);
text(locs1+.02,pks1,num2str((1:numel(pks1))'))

figure(1)
subplot(413)
baseTime = 0;

% beat times in seconds
beatTimes = baseTime + (1/fs)*(locs'-1);

%attack_time = round(locs);
bufferLen = size(x,1);
sampleLen = size(y,1);

sampleTrack = zeros(bufferLen+sampleLen,2);

% Place sample at each onset
for k = 1:numel(locs)
   %if index is < 0, i.e. -20, then cut off 20 samples, but also cut off 
   %20 samples from y
   %idx+sampleLen-1 > length(sampleTrack)
   
   difflocs = locs(k)/locs1;
   
   if(difflocs > 0)
      idx = locs(k) - difflocs;
   else
      idx = locs(k);
   end
   
   %if(idx < 0)
       %sampleTrack(idx : idx+sampleLen-1, :) = sampleTrack(find(sampleTrack>0,1):end);   
   %end
   
   sampleTrack(idx : idx+sampleLen-1, :) = sampleTrack(idx : idx+sampleLen-1, :) + y;
   
end

%a_short = sampleTrack(find(sampleTrack>0,1):end);

Guide = x(44100*1 : 44100*7);
Aligned = sampleTrack(44100*1 : 44100*7);

mix = Guide + Aligned;

ax = 0:N-overlap:(N-overlap) * (length(HFC) - 1);
plot(sampleTrack(:,1))
colorbar
xlabel('Trigger Placement')
axis([min(ax),max(ax),min(sampleTrack(:,1)),max(sampleTrack(:,1))])

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
xlabel('Spectral Flux')
axis([0,max(ax),0,max(spectral_flux)])

[pks2,locs2] = findpeaks(spectral_flux,ax,'MinPeakDistance',.1,'MinPeakHeight',75);
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
ax = 0:N-overlap:(N-overlap) * (length(HFC) - 1);
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
ax=linspace(0,length(z)/fs,length(z));
plot(ax,z);
xlabel('Mistrigger')

figure(2)
subplot(313)
ax=linspace(0,length(y)*fs,length(y));
plot(ax,y)
xlabel('Clean Sample')