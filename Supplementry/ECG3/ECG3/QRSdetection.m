function [ Rcoor] = QRSdetection( ECG,fs)
%QRS localization OUTPUT: Rposition
N = length(ECG);
tm = 1/fs:1/fs:N/fs;
x = 4;
y = 5;

wt = modwt(ECG,'db6',y);
wtrec = zeros(size(wt));
wtrec(x:y,:) = wt(x:y,:);
y = imodwt(wtrec,'db6');
y = abs(y).^2;
th = 0.15*max(y);
[qrspeaks,locs] = findpeaks(y,tm,'MinPeakHeight',th,...
    'MinPeakDistance',0.150);
Rcoor = round(locs*fs);
%%%%% the following part correction for faulse detection %%%%%%
windowQRS = [ceil(fs*0.05) ceil(fs*0.1)];
for i = 1:length(Rcoor)  
    if ECG(Rcoor(i)) <0
        % exclude the incomlete cardiac cycle
        QRSon = max( [Rcoor(i) - windowQRS(1) 1]);
        QRSoff = min( [Rcoor(i) + windowQRS(2) length(ECG)] );
        R = ECG(QRSon:QRSoff);
        [M,I] = max(R); % R peak is supposed to be the maxima within QRS complex
        Rcoor(i) = min( [ I+QRSon length(ECG)]);
    end
end

% figure;plot(tm,ECG);hold on;scatter(locs,ECG(Rcoor),'ro')
% xlabel('Seconds'); ylabel('Amplitude')
% title('ECG R peaks detection')
end

