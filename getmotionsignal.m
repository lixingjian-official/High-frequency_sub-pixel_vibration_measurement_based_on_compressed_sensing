function [S] = getmotionsignal(dataDir,index, varargin)
nscales = 3;
norientations = 2;
if mod(norientations,2) == 1
    error('The number of norientations can not be odd !!!')
end

tic;
startTime = toc;
%% Parameters
defaultnframes = 0;
defaultDownsampleFactor = 1;
defaultsamplingrate = -1;
p = inputParser();
addOptional(p, 'NFrames', defaultnframes, @isnumeric);  
addOptional(p, 'DownsampleFactor', defaultDownsampleFactor, @isnumeric);   
addOptional(p, 'SamplingRate', defaultsamplingrate, @isnumeric);   
parse(p, varargin{:});
nScales = nscales;
nOrients = norientations;
dSampleFactor = p.Results.DownsampleFactor;
samplingrate = p.Results.SamplingRate;


'Reading reference frame'
refFrame = imread([dataDir, '\','refFrame.tif']);
'Successfully read reference frame'

if(dSampleFactor~=1)
    refFrame = imresize(refFrame,dSampleFactor);
end

nF = length(index);

%% Build Complex-Valued Steerable Pyramid
[pyrRef, pind] = buildSCFpyr(refFrame, nScales, nOrients-1); %Complex-Valued Steerable Pyramid
signalffs = zeros(nScales,nOrients,nF);%All Frame Scale and Orients 
ampsigs = zeros(nScales,nOrients,nF);

%% Process
nF

for q = 1:nF
    if(mod(q,floor(nF/100))==1) %
        progress = q/nF;
        currentTime = toc;
        ['Progress:' num2str(progress*100) '% done after ' num2str(currentTime-startTime) ' seconds.']
    end
    vframein = imread([dataDir, '\',num2str(index(q)),'.tif']);
    if(dSampleFactor == 1)
        fullFrame = vframein;
    else
        fullFrame = imresize(vframein,dSampleFactor);
    end            
    im = fullFrame;
    pyr = buildSCFpyr(im, nScales, nOrients-1);
    pyrAmp = abs(pyr);
    pyrDeltaPhase = mod(pi+angle(pyr)-angle(pyrRef), 2*pi) - pi; 
    for j = 1:nScales    
        for k = 1:nOrients
            bandIdx = 1 + (j-1)*nOrients + k;
            amp = pyrBand(pyrAmp, pind, bandIdx); 
            phase = pyrBand(pyrDeltaPhase, pind, bandIdx); 
            phasew = phase.*(abs(amp).^2);      
            sumamp = sum(abs(amp(:)));
            ampsigs(j,k,q)= sumamp;
            signalffs(j,k,q)=mean(phasew(:))/sumamp; 
        end
    end    
end

S.samplingRate = samplingrate;

%% Vibration Signals in the X Direction and Y Direction
sigOut_horizontal = zeros(nF, 1);
sigOut_vertical = zeros(nF, 1);
for q=1:nScales
    for p=1:nOrients
        sigaligned = squeeze(signalffs(q,p,:));
        if mod(p,2) == 1
            sigOut_horizontal = sigOut_horizontal-sigaligned;
        else
            sigOut_vertical = sigOut_vertical+sigaligned;
        end

    end
end
S.horizontal = sigOut_horizontal;
S.vertical = sigOut_vertical;
S.average = mean(reshape(double(signalffs),nScales*nOrients,nF)).';
end
