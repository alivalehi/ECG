clear
load('nMITDS1_4.mat')

N = [nMITDS1_4(find(nMITDS1_4(:,5)==1),:) find(nMITDS1_4(:,5)==1)]; 
V = [nMITDS1_4(find(nMITDS1_4(:,5)==2),:) find(nMITDS1_4(:,5)==2)];
S = [nMITDS1_4(find(nMITDS1_4(:,5)==3),:) find(nMITDS1_4(:,5)==3)];
F = [nMITDS1_4(find(nMITDS1_4(:,5)==4),:) find(nMITDS1_4(:,5)==4)];

figure;lscatter(N(1:50:end,1),N(1:50:end,2),N(1:50:end,6),'TextColor','b');
hold on;lscatter(V(1:10:end,1),V(1:10:end,2),V(1:10:end,6),'TextColor','r');
hold on;lscatter(S(1:2:end,1),S(1:2:end,2),S(1:2:end,6),'TextColor','m');
hold on;lscatter(F(1:end,1),F(1:end,2),F(1:end,6),'TextColor','c');
title('prejection on 2-D by PCA')
%%
figure;
namelist1=[101;106;108;109;112;114;115;116;118;119;122;124;201;203;205;207;208;209;215;220;223;230];
counter = 1;
searching = [12199 11077 5583 10687];
lab = ['N' 'V' 'S' 'F'];
for s = 1:4
    counter = 1;
for i = 1:length(namelist1)
    testName = ['\\EGRSHARES\Homes\NAU-STUDENTS\jc3464\Documents\MATLAB\ECG2\MITDB(non-normalized)\' num2str(namelist1(i)) 'Interval.mat'];
    load(testName);
    if counter<searching(s)&& searching(s)<(counter + size(Ind,1)-1)
        disp(['data point found in patient #',num2str(namelist1(i)), '"s record'])
        start = 1;
        while (start+counter-1)~=searching(s)
            start = start+1;
        end
        found = Ind(start,:);
        break
    end
    counter = counter+size(Ind,1)-1;
end

matName = ['\\EGRSHARES\Homes\NAU-STUDENTS\jc3464\Documents\MATLAB\ECG\online database\MITDBmat\' num2str(namelist1(i)) 'm.mat'];
ECGsignal = load(matName);
%figure;
if lab(s) == 'N'
    subplot(2,2,1);
    plot(ECGsignal.val(found(1):found(2)));title(['patient#',num2str(namelist1(i)),' label:',lab(s),'data number',num2str(searching(s))])
elseif lab(s) == 'V'
    subplot(2,2,2);
    plot(ECGsignal.val(found(1):found(2)));title(['patient#',num2str(namelist1(i)),' label:',lab(s),'data number',num2str(searching(s))])
elseif lab(s) == 'S'
    subplot(2,2,3);
    plot(ECGsignal.val(found(1):found(2)));title(['patient#',num2str(namelist1(i)),' label:',lab(s),'data number',num2str(searching(s))])
else
    subplot(2,2,4);
    plot(ECGsignal.val(found(1):found(2)));title(['patient#',num2str(namelist1(i)),' label:',lab(s),'data number',num2str(searching(s))])
end
end
