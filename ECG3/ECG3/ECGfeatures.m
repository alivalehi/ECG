function [ Features,varagout ] = ECGfeatures( ECGsignal,fs,baseline,gain,window,step)
%This is the function which combines all the feature extraction fct for ECG
%signal.
%   Input: raw ECG signal, Smapling frequency (906Hz)
%   Output: the features list for further machine learning part
%           Pposition,Qposition,Rposition,Sposition,Tposition, median of
%           RR interval, QRS duration, 

%% Pre processing, denoise,smoothing
ECG = ECGsignal;%(5:N-5);
clear ECGsignal
[ ECG_DENO ] = Preprocessing( ECG,fs,baseline,gain);
%% Get the position of PQRST pour cet echantillon

[ Rcoor ] = QRSdetection( ECG_DENO,fs);
[ECGpeaks] = QSpeaks( ECG_DENO,Rcoor,fs );
% figure;plot(ECG);
% hold on;
% for i = 1:length(Rposition)
%     R(i) = ECG(Rposition(i));
% end
% plot(Rposition,R,'*');
[ Features,ind ] = SegmentFeature2( ECG_DENO,ECGpeaks,fs,window,step );
varagout = ind;
%[Pposition,QRS_ON, Qposition, Sposition ,QRS_OFF, Tposition] = QSpeaks( ECG_DENO,Rposition,fs );
%[ Features ] = SegmentFeature(ECG_DENO, Rposition,Pposition, Qposition, Sposition , Tposition,QRS_ON,QRS_OFF,fs );
%csvwrite('MITstdb1min.csv', Features);
%NormFeatures = zscore(Features);

end

