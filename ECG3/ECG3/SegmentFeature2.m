function [ Features,varagout ] = SegmentFeature2( ECG_DENO,ECGpeaks,fs, varargin )
% Second Version of Segment Features
%   Initialization
WindowWidth = varargin{1};
StepLength = varargin{2};%  = 3;
Rposition = ECGpeaks(:,4);
IndValidSeg = [];
Segment = [];
Features = [];
h = length(ECG_DENO)/length(Rposition);
for Rcentral = ceil(WindowWidth/2):StepLength:length(Rposition)-floor(WindowWidth/2)
    Rs = Rcentral-floor(WindowWidth/2);
    Re = Rcentral+floor(WindowWidth/2);
    segmentPKS = ECGpeaks(Rs:Re,:);
    Rstart = Rposition(Rs);
    Rend = Rposition(Re);
    RRthisCycle = [];
    for j = Rs:Re-1
        RRthisCycle(end+1) = Rposition(j+1)-Rposition(j);
    end
    %rr = RRthisCycle./fs;
    if  Rstart - round(h/3) >= 1 && Rend + round(h*2/3) <= length(ECG_DENO)
        IndValidSeg(end+1) = segmentPKS(1,1);
        IndValidSeg(end+1) = segmentPKS(WindowWidth,7);
        %the index boundary of this segment
         Segment = ECG_DENO(Rstart - round(h/3):Rend + round(h*2/3));
         thisPKS = segmentPKS - (Rstart - round(h/3)) + 1;
         %to make sure the peak indexes won't exceed the range
         if thisPKS(1,1) <1
             thisPKS(1,1) = 1;
         end
         if thisPKS(WindowWidth,7) > length(Segment)
             thisPKS(WindowWidth,7) = length(Segment);
         end         
        TimePerSeg = zeros(3,WindowWidth-1);
        FreqPerSeg = zeros(4,WindowWidth-1);
        WavePerSeg = zeros(3,WindowWidth-1);
        MorphPerSeg = zeros(1,WindowWidth-1);
        for indCycle = 1:WindowWidth
            if indCycle == 1
                Cend = max([thisPKS(indCycle,7) thisPKS(1,4)+ round(2*(thisPKS(2,4)-thisPKS(1,4))/3)]);
                thisCycle = Segment(1:Cend);
                peaks = thisPKS(indCycle,:);
            elseif indCycle == WindowWidth
                Cbegin =  min([thisPKS(indCycle,1) thisPKS(WindowWidth,4)- round((thisPKS(WindowWidth,4)-thisPKS(WindowWidth-1,4))/3)]);
                thisCycle = Segment(Cbegin:length(Segment));
                peaks = thisPKS(indCycle,:) - Cbegin+1;
            else
                Cbegin = min([thisPKS(indCycle,1) thisPKS(indCycle,4)- round((thisPKS(indCycle,4)-thisPKS(indCycle-1,4))/3)]);
                Cend = max([thisPKS(indCycle,7) thisPKS(indCycle,4)+ round(2*(thisPKS(indCycle+1,4)-thisPKS(indCycle,4))/3)]);
                thisCycle = Segment(Cbegin:Cend);
                peaks = thisPKS(indCycle,:) - Cbegin+1;
            end
            if all(peaks>0)
                deb =0;
            else
                deb = 1;
            end
            %correction
            if peaks(1) < 1
                peaks(1) = 1;
            elseif peaks(7) > length(thisCycle)
                peaks(7) = length(thisCycle);
            end
            [ Time,Freq,Wave,Morph ] = CycleFeature( thisCycle,peaks(1),peaks(2),peaks(3),peaks(4),peaks(5),peaks(6),peaks(7),fs);
             TimePerSeg(:,indCycle) = Time;
             FreqPerSeg(:,indCycle) = Freq;
             WavePerSeg(:,indCycle) = Wave;
             MorphPerSeg(:,indCycle) = Morph;
        end
        [ TimeSta,FreqSta,WaveSta,MorphSta ] = SegmentStatistics( TimePerSeg, FreqPerSeg,WavePerSeg,MorphPerSeg ); 
        SAE = sum(Segment.^2)*fs/length(Segment);
        MPP = max(Segment);
        MNP = -min(Segment);
        PER = max(Segment)*length(Segment)/(fs*sum(Segment.^2));
        Feature = [reshape(TimeSta,1,[]) reshape(FreqSta,1,[]) reshape(WaveSta,1,[]) reshape(MorphSta,1,[]) SAE MPP MNP PER];
        Feature = [mean(RRthisCycle)/fs (mean(RRthisCycle-h))/fs Feature];
        Features = [Features;Feature];       
    end
        
end
% X = ['we have ',num2str(length(IndValidSeg)),' valid segment'];
% disp(X)
IndSeg = zeros(2,length(IndValidSeg)/2);
for ii = 1:length(IndValidSeg)
    if mod(ii,2) == 1
        IndSeg(1,ceil(ii/2)) = IndValidSeg(ii);
    else
        IndSeg(2,ii/2) = IndValidSeg(ii);
    end
end
varagout = IndSeg';
end

