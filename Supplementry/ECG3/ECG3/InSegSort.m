function [ PerCycle ] = InSegSort( Segment,P,QRS_ON, Q,R, S ,QRS_OFF, T )
%UNTITLED34 Summary of this function goes here
%   Detailed explanation goes here
L = [length(P),length(QRS_ON),length(Q),length(R),length(S),length(QRS_OFF),length(T)];
N = max(L);
if range(L) >=2
    PerCycle = 0;
    return
end
if P(1)>QRS_ON(1)
    [P1V,P1] = max(Segment(1:QRS_ON(1)));
    P = [P1 P];
end
if T(end)<QRS_OFF(end)
    [TV,Te] = max(Segment(QRS_OFF(end):end));
    T = [T (Te+QRS_OFF(end))];
end
PerCycle = zeros(N,length(L));
for i = 1:N
    PerCycle(i,:) = [P(i) QRS_ON(i) Q(i) R(i) S(i) QRS_OFF(i) T(i)];
    if issorted(PerCycle(i,:)) ==0
        disp('abnormal heartbeat')
        PerCycle(i,:) = sort(PerCycle(i,:));
    end
end
    

end

