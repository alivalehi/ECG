ds = input('DS#: ');
if ds == 1
clear;
namelist1=[101;106;108;109;112;114;115;116;118;119;122;124;201;203;205;207;208;209;215;220;223;230];
%patient_record = zeros(3700,22);
%filename = 'patient record.xlsx';
%load('MITDS1.mat')
MITDS1 = [];
Vdb = [];
Sdb = [];
Fdb = [];
for o = 1:length(namelist1)
%for o = 6:7
eval(['annotName = ''\\EGRSHARES\Homes\NAU\jc3464\Documents\MATLAB\ECG\online database\MITDBannot\' num2str(namelist1(o)) '.xlsx'';'])
[num,txt,raw] = xlsread(annotName);
colHeadings = {'Time' 'Sample' 'Type'};
readAnno = cell2struct(raw(2:end,1:3), colHeadings, 2);
beat = size(raw,1)-1;
Label = zeros(beat,2);
for r = 1:beat
    if strcmp(readAnno(r).Type,'N') || strcmp(readAnno(r).Type,'e') || strcmp(readAnno(r).Type,'j') || strcmp(readAnno(r).Type,'L') ||strcmp(readAnno(r).Type,'R')
        Label(r,1) = readAnno(r).Sample;
        Label(r,2) = 1;
    elseif strcmp(readAnno(r).Type,'V') || strcmp(readAnno(r).Type,'E')
        Label(r,1) = readAnno(r).Sample;
        Label(r,2) = 2;
    elseif strcmp(readAnno(r).Type,'a') || strcmp(readAnno(r).Type,'A') || strcmp(readAnno(r).Type,'J') || strcmp(readAnno(r).Type,'S')
        Label(r,1) = readAnno(r).Sample;
        Label(r,2) = 3;
    elseif strcmp(readAnno(r).Type,'F')
        Label(r,1) = readAnno(r).Sample;
        Label(r,2) = 4;
    end
end
Lab = Label(find(Label(:,2)~=0),:);
matName = ['\\EGRSHARES\Homes\NAU\jc3464\Documents\MATLAB\ECG\online database\MITDBmat\' num2str(namelist1(o)) 'm.mat'];
ECGsignal = load(matName); fs = 360;
ECGsignal = ECGsignal.val;
baseline = 0;gain = 1024;
%patient_record(1:length(Lab),o) = Label(Lab,2);
%xlswrite(filename,patient_record);
wind = 1;
stp = 1;
[ Features,Ind ] = ECGfeatures( ECGsignal,fs,baseline,gain,wind,stp );
featureFile = [num2str(namelist1(o)) 'NF.mat']; intervalFile = [num2str(namelist1(o)) 'Interval.mat'];
labelFile = [num2str(namelist1(o)) 'Label'];
snum = size(Features,1);
Features = [Features zeros(snum,1)];
CycleLab = zeros(size(Ind,1),wind);
j = 1;
while j <= size(Ind,1)
    for i = 1:size(Lab,1)-wind+1
        if Lab(i,1)> Ind(j,1) && Lab(i+wind-1,1)<Ind(j,2)            
            CycleLab(j,:) = Lab(i:i+wind-1,2);
            nvsf = [1 2 3 4];
            res = double(ismember(nvsf,(CycleLab(j,:))));
            if isequal(res,[1 0 0 0])
               Features(j,29) = 1;
            elseif isequal(res,[1 1 0 0]) || isequal(res,[0 1 0 0])
                Features(j,29) = 2;
                Vdb = [Vdb;Features(j,:)];
            elseif isequal(res,[1 0 1 0]) || isequal(res,[0 0 1 0])
                Features(j,29) = 3;
                Sdb = [Sdb;Features(j,:)];
            elseif isequal(res,[1 0 0 1]) || isequal(res,[0 1 0 1]) || isequal(res,[0 0 0 1])
                Features(j,29) = 4;
                Fdb = [Fdb;Features(j,:)];
            end
            
        end        
    end
    j =j+1;
end
cc = find(Features(:,29)~=0);
Features = Features(cc,:);
save(featureFile,'Features');save(intervalFile,'Ind');save(labelFile,'Label');
MITDS1 = [MITDS1; Features];
end
nMITDS1 = zscore(MITDS1);
save('mitds1','MITDS1')
save('nmitds1','nMITDS1')
save('mitds1v','Vdb')
save('mitds1s','Sdb')
save('mitds1f','Fdb')

%%
else
clear
MITDS2 = [];
namelist2=[100;103;105;111;113;117;121;123;200;202;210;212;213;214;219;221;222;228;231;232;233;234];
for o = 1:length(namelist2)
eval(['annotName = ''\\EGRSHARES\Homes\NAU-STUDENTS\jc3464\Documents\MATLAB\ECG\online database\MITDBannot\' num2str(namelist2(o)) '.xlsx'';'])
[num,txt,raw] = xlsread(annotName);
colHeadings = {'Time' 'Sample' 'Type'};
readAnno = cell2struct(raw(2:end,1:3), colHeadings, 2);
beat = size(raw,1)-1;
Label = zeros(beat,2);
for r = 1:beat
    if strcmp(readAnno(r).Type,'N') || strcmp(readAnno(r).Type,'e') || strcmp(readAnno(r).Type,'j') || strcmp(readAnno(r).Type,'L') ||strcmp(readAnno(r).Type,'R')
        Label(r,1) = readAnno(r).Sample;
        Label(r,2) = 1;
    elseif strcmp(readAnno(r).Type,'V') || strcmp(readAnno(r).Type,'E')
        Label(r,1) = readAnno(r).Sample;
        Label(r,2) = 2;
    elseif strcmp(readAnno(r).Type,'a') || strcmp(readAnno(r).Type,'A') || strcmp(readAnno(r).Type,'J') || strcmp(readAnno(r).Type,'S')
        Label(r,1) = readAnno(r).Sample;
        Label(r,2) = 3;
    elseif strcmp(readAnno(r).Type,'F')
        Label(r,1) = readAnno(r).Sample;
        Label(r,2) = 4;
    end
end
Lab = Label(find(Label(:,2)~=0),:);
matName = ['\\EGRSHARES\Homes\NAU-STUDENTS\jc3464\Documents\MATLAB\ECG\online database\MITDBmat\' num2str(namelist2(o)) 'm.mat'];
ECGsignal = load(matName); fs = 360;
ECGsignal = ECGsignal.val;
baseline = 0;gain = 1024;
%patient_record(1:length(Lab),o) = Label(Lab,2);
%xlswrite(filename,patient_record);
wind = 1;
stp = 1;
[ Features,Ind ] = ECGfeatures( ECGsignal,fs,baseline,gain,wind,stp );
featureFile = [num2str(namelist2(o)) 'NF.mat']; intervalFile = [num2str(namelist2(o)) 'Interval.mat'];
labelFile = [num2str(namelist2(o)) 'Label'];
snum = size(Features,1);
Features = [Features zeros(snum,1)];
CycleLab = zeros(size(Ind,1),wind);
j = 1;
while j <= size(Ind,1)
    for i = 1:size(Lab,1)-wind+1
        if Lab(i,1)> Ind(j,1) && Lab(i+wind-1,1)<Ind(j,2)            
            CycleLab(j,:) = Lab(i:i+wind-1,2);
            nvsf = [1 2 3 4];
            res = double(ismember(nvsf,(CycleLab(j,:))));
            if isequal(res,[1 0 0 0])
               Features(j,29) = 1;
            elseif isequal(res,[1 1 0 0]) || isequal(res,[0 1 0 0])
                Features(j,29) = 2;
            elseif isequal(res,[1 0 1 0]) || isequal(res,[0 0 1 0])
                Features(j,29) = 3;
            elseif isequal(res,[1 0 0 1]) || isequal(res,[0 1 0 1]) || isequal(res,[0 0 0 1])
                Features(j,29) = 4;
            end            
        end        
    end
    j =j+1;
end
cc = find(Features(:,29)~=0);
Features = Features(cc,:);
save(featureFile,'Features');save(intervalFile,'Ind');save(labelFile,'Label');
MITDS2 = [MITDS2; Features];
end
save('mitds2','MITDS2')
end