addpath /Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/util
%Sound files to be decomposed
mixDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/mixes/');
mixFileList = getFileNames(mixDirectory ,'wav');
%Directories of isolated snare sound files%
snareDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/close mics/snares/');
snareFileList = getFileNames(snareDirectory ,'wav');
%Directories of isolated kick sound files%
kickDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/close mics/kicks/');
kickFileList = getFileNames(kickDirectory ,'wav');
%Directories of isolated high hat sound files%
hhDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/close mics/hh/');
hhFileList = getFileNames(hhDirectory ,'wav');
%Mat file directory%
closematFileDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/closematFiles/');
closematFileList = getFileNames(closematFileDirectory ,'mat');
%Mat file directory, mixes%
mixmatFileDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/mixmatFiles/');
mixmatFileList = getFileNames(mixmatFileDirectory ,'mat');
%TestMat file directory, mixes%
testmixmatFileDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/testmixmatFiles/');
testmixmatFileList = getFileNames(testmixmatFileDirectory ,'mat');
%TestMat file directory, close mics%
testclosematFileDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/testclosematFiles/');
testclosematFileList = getFileNames(testclosematFileDirectory ,'mat');

%% NMFD Decomp and save on mix files %%
for i=1:length(mixFileList)
    %% parameters %%
    %length of the FFT%
    Nfft = 512;
    %number of time-frequency atoms%
    R = 3;
    % time-length of the atoms%
    T = 120;
    % number of iterations%
    Niter = 20;
    %% preparing the data %%
    disp('reading mix file...')
    %Load the sound file to be decomposed%
    [s,sr] = audioread([mixDirectory mixFileList{i}]);
    %convert to mono by default%
    s = toMono(s);
    %compute the spectrograms of mix file and close mics%
    sp = stft(s,Nfft,hamming(Nfft,'periodic'),Nfft/4);
    V = abs(sp);
    M = size(V,1);
    N = size(V,2);
    %first decomposition%
    [W,H] = NMFD(V,R,T,Niter);
    %post processing%
    initVal.W = W;
    initVal.H = max(H,max(H(:))/10)-max(H(:))/10+0.00001;
    %second decomposition%
    [W,H] = NMFD(V,R,T,10,initVal);
    %separation of the sounds of each component%
    Lambda = cell(R,1);
    for z = 1:R
        Lambda{z} = zeros(M,N);
        for f = 1:M
            v = reshape(W(f,z,:),T,1);
            cv = conv(v,H(z,:));
            Lambda{z}(f,:) = Lambda{z}(f,:) + cv(1:N);
        end
    end

    LambdaTot = zeros(M,N);
    for z = 1:R
        LambdaTot = LambdaTot +Lambda{z};
    end

    for z = 1:R
        xs{z} = istft(sp.*Lambda{z}./LambdaTot,Nfft,hamming(Nfft,'periodic')',Nfft/4);
        %soundsc(xs{z},sr);
    end

    filename = sprintf('%s_%d','nmfdWH',i);

    disp('storing mix file... ')
    save([mixmatFileDirectory filename '.mat'],'W','H','Lambda','xs')

end

%% Saving Close Mic Files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(snareFileList)
    
    disp('reading close files...')
    
    [snclose,fs]=audioread([snareDirectory snareFileList{i}]); 
    [kiclose,fs]=audioread([kickDirectory kickFileList{i}]);  
    [hhclose,fs]=audioread([hhDirectory hhFileList{i}]);
    
    disp('extracting from close files...')
    
    snclose = toMono(snclose); % mono
    kiclose = toMono(kiclose); % mono
    hhclose = toMono(hhclose); % mono
    
    filename = sprintf('%s_%d','closemics',i);
    
    disp('storing close files... ')
    save([closematFileDirectory filename '.mat'],'snclose','kiclose','hhclose')
    
end

% Pull out the mat files to the test files
%Set 1: 5, 6, 7, 28, 36, 38, 41, 43, 49
%       53, 54, 56, 58, 59, 70, 75
%       79, 84, 87, 93, 97
%20 numbers randomly selected between 1 and 100

%% Reload them into Matlab / SVM

% Main Files
[features,labels,HofLabel] = LoadingProcedure(mixmatFileList,mixmatFileDirectory,closematFileList,closematFileDirectory);
% Test Files
[testfeatures,testlabels,testHofLabel] = LoadingProcedure(testmixmatFileList,testmixmatFileDirectory,testclosematFileList,testclosematFileDirectory);


% save H of labels (locations of 'notes' per label)
disp('storing H of label... ')
save('H_of_label.mat','HofLabel')
% save test H of labels (locations of 'notes' per label)
disp('storing test H of label... ')
save('test_H_of_label.mat','testHofLabel')

% save features and labels
disp('storing features and labels... ')
save('features_labels.mat','features','labels')
% save features and labels
disp('storing features and labels... ')
save('test_features_labels.mat','testfeatures','testlabels')

% load features and labels, put them in model
load features_labels
X = features;
Y = categorical(labels);
classOrder = unique(Y);
rng(1); % For reproducibility

% do this to the 80 left after pulling 20
t = templateSVM('Standardize',1);
CVMdl = fitcecoc(X,Y,'Learners',t,'ClassNames',classOrder);
CMdl = CVMdl;

% Extract trained, compact classifier
%testInds = test(CVMdl.Partition);  % Extract the test indices
%save file names as it goes
%save('testInds.mat','testInds');

% fileNum = ceil(find(testInds) / 3);
% componentNum = mod(find(testInds)-1, 3) + 1;

% Get XTest/YTest from reloading the test files
% Indices from 20 pulled out
% do loading procedure again for test files. Make it a function

load test_features_labels
% XTest = features, YTest = Labels
XTest = testfeatures;
YTest = testlabels';

preds = predict(CMdl,XTest);
table(YTest,preds,...
    'VariableNames',{'TrueLabels','PredictedLabels'})


%Read it from the saved spread sheet%
TrueLabels = readtable('TestInds_Number_Name_Final.xlsx','Range','E1:E61');
PredictedLabels = readtable('TestInds_Number_Name_Final.xlsx','Range','F1:F61');
%convert table to variable array for operations%
TrueLabels = table2array(TrueLabels);
PredictedLabels = table2array(PredictedLabels);
%save new variables arrays%
save('true_predicted_labels','TrueLabels','PredictedLabels');
%create confusion matrix%
g1 = TrueLabels';		% Known groups
g2 = PredictedLabels';	% Predicted groups
[C,order] = confusionmat(g1,g2);
%Look at confusion matrix results
C
order

% Across the top row is kick, 17 correct, though 4 were snares, thought 3
% were high hats i.e.
% 










%%%%%%%%%%%%%%%%%%%%%%MIDI Section%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
d = length(s)/sr;
nums = (1:N)/N*d;
[pksH1,locsH1] = findpeaks(H(1,:),nums,'MinPeakDistance',.1,'MinPeakHeight',.1);
[pksH2,locsH2] = findpeaks(H(2,:),nums,'MinPeakDistance',.1,'MinPeakHeight',.1);
[pksH3,locsH3] = findpeaks(H(3,:),nums,'MinPeakDistance',.1,'MinPeakHeight',.1);

HFCKick=zeros(1,size(spkclose,2));
for k = 2:size(spkclose,2)
    for j = 2 : size(spkclose,1)
        HFCKick(k) = HFCKick(k)+(abs(spkclose(j,k))).^2 * j;
    end
end

HFCSnare=zeros(1,size(spsclose,2));
for k = 2:size(spsclose,2)
    for j = 2 : size(spsclose,1)
        HFCSnare(k) = HFCSnare(k)+(abs(spsclose(j,k))).^2 * j;
    end
end

HFCHH=zeros(1,size(sphclose,2));
for k = 2:size(sphclose,2)
    for j = 2 : size(sphclose,1)
        HFCHH(k) = HFCHH(k)+(abs(sphclose(j,k))).^2 * j;
    end
end

[pkskiclose,locskiclose] = findpeaks(HFCKick,nums,'MinPeakDistance',.5,'MinPeakHeight',200);
[pkssnclose,locssnclose] = findpeaks(HFCSnare,nums,'MinPeakDistance',.1,'MinPeakHeight',1000);
[pkshhclose,locshhclose] = findpeaks(HFCHH,nums,'MinPeakDistance',.3,'MinPeakHeight',2000);
    










%% plot results

d = length(s)/sr;

figure 
imagesc((1:N)/N*d,0:0.05:sr/2000,db(V))
title('Original Spectrogram')
xlabel('time (s)')
ylabel('frequency (kHz)')

axis xy
mx = max(caxis);
caxis([mx-120,mx])
dyn = caxis;

figure
subplot(131)
imagesc((1:T)/N*d,0:0.05:sr/2000,db(reshape(W(:,1,:),M,T)))
caxis(dyn)
axis xy
title('First template')
xlabel('time (s)')
ylabel('frequency (kHz)')

subplot(132)
imagesc((1:T)/N*d,0:0.05:sr/2000,db(reshape(W(:,2,:),M,T)))
caxis(dyn)
axis xy
title('Second template')
xlabel('time (s)')
ylabel('frequency (kHz)')

subplot(133)
imagesc((1:T)/N*d,0:0.05:sr/2000,db(reshape(W(:,3,:),M,T)))
caxis(dyn)
axis xy
title('Third template')
xlabel('time (s)')
ylabel('frequency (kHz)')

figure
subplot(311)
plot((1:N)/N*d,(H(1,:)))
title('Activation of the 1st template')
xlabel('time (s)')


subplot(312)
plot((1:N)/N*d,(H(2,:)))
title('Activation of the 2nd template')
xlabel('time (s)')
subplot(313)
plot((1:N)/N*d,(H(3,:)))
title('Activation of the 3rd template')
xlabel('time (s)')