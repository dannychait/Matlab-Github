function [features,labels,HofLabel] = LoadingProcedure(mixmatFileList,mixmatFileDirectory,closematFileList,closematFileDirectory)

Nfft = 512;
R=3;

j=1;
for i=1:length(mixmatFileList)
    
    load([mixmatFileDirectory mixmatFileList{i}])
    %savedNmfFiles = load(matFilename);
    %savedNmfFiles = struct2cell(savedNmfFiles);
    
    load([closematFileDirectory closematFileList{i}])
    
    %matclosemicFilename=strrep(closemicFilename,'.wav','.mat');
    %save(matclosemicFilename,'kiclose','snclose','hhclose')
    %closeMicFiles = load(matclosemicFilename);
    %closeMicFiles = struct2cell(closeMicFiles);
    
    spsclose = stft(snclose,Nfft,hamming(Nfft,'periodic'),Nfft/4);
    spkclose = stft(kiclose,Nfft,hamming(Nfft,'periodic'),Nfft/4);
    sphclose = stft(hhclose,Nfft,hamming(Nfft,'periodic'),Nfft/4);
 
    closeMicNames = {'kick','snare','high hat'};
    matches = {'kick','snare','high hat'}; %Default Matching

%Match Spectrogram Function  
    %use min function on l
    %Once you know the lambda{z} the drum is closest to
    %put the string name of the drum in matches{z}
    %z in Lambda is the same as z in matches
    
    %Try different values of R
    %Try different songs
    %Try HFC
    %griffin lim
    
    l = zeros(R,R);
    for z = 1:R
        l(1,z) = sum(sum((db(spkclose) - db(Lambda{z})).^2,1),2);
        l(2,z) = sum(sum((db(spsclose) - db(Lambda{z})).^2,1),2);
        l(3,z) = sum(sum((db(sphclose) - db(Lambda{z})).^2,1),2);
    end
 
    for z = 1:R
        [~,ind] = min(l(:,z));
        matches{z} = closeMicNames{ind};  
    end
   
    for f = 1:R
        features(j,:) = reshape(W(:,:,f), 1, []);
        %save H peaks at same index
        HofLabel{j} = H(f,:);
        
        labels{j} = matches{f};
        
     
        j = j + 1;
   
 
    end
end
