%%%% Sytem of classification May 1st 2018
%%%% transformation to be implemented, transformation before classification

% loading data base(features)
test = 1; % validate or test
DS1 = [101;106;108;109;112;114;115;116;118;119;122;124;201;203;205;207;208;209;215;220;223;230];
DS2 = [100;103;105;111;113;117;121;123;200;202;210;212;213;214;219;221;222;228;231;233;234];
if test == 1
    record_list = DS2;
else
    record_list = DS1;
end

rec_num = length(record_list);

for r = 1:rec_num
% load individual record
filename = ['\\EGRSHARES\Homes\NAU\jc3464\Documents\MATLAB\ECG2\MITDB(non-normalized)\' num2str(record_list(r)) 'NF.mat'];
data = load(filename);
DD = size(data,2)-1; % dimension of data (should be 28)
NN = size(data,1);% length of data
%forming initial N cluster (first b20% of the total N samples)
N_all = find(data(:,DD+1)==1);


end