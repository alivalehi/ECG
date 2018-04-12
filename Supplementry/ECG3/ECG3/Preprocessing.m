function ECG = Preprocessing(ECG_RAW,fs,baseline,gain)
fmit = 360;
%%% baseline removal and normalize
ECG_RAW = ECG_RAW - baseline;
ECG_RAW = ECG_RAW./gain;
opol = 5;
t = (1:length(ECG_RAW));
[p,s,mu] = polyfit(t,ECG_RAW,opol);
f_y = polyval(p,t,[],mu);
ecg = ECG_RAW - f_y;

%%% wavelet denoising
lev = 1+round(log2(fs/fmit));
ecg = wden(ecg,'heursure','s','mln',lev,'db8'); % \cite{optidb8}

%%% resample to 360Hz
ECG = resample(ecg,fmit,fs);
end

