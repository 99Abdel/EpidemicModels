function Times = timeChanges(t1,t2,I)

tmax = t2-t1+1;
dt = 1;

t = 1:dt:tmax;
[asseX, Data] = smoothCurve(t,I);

Distance = 50;
Prominence = 500;

[Maxima,MaxIdx] = findpeaks(Data,asseX,'MinPeakDistance',Distance,'MinPeakProminence',Prominence);
DataInv = 1.01*max(Data) - Data;
[Minima,MinIdx] = findpeaks(DataInv,asseX,'MinPeakDistance',Distance,'MinPeakProminence',Prominence);
Minima = Data(MinIdx);


Times = sort([MinIdx MaxIdx]);

end


% figure
% findpeaks(Data,asseX,'MinPeakDistance',Distance,'MinPeakProminence',Prominence);
% text(MaxIdx+.02,Maxima,num2str((1:numel(Maxima))'))
% figure
% findpeaks(DataInv,asseX,'MinPeakDistance',Distance,'MinPeakProminence',Prominence);
% text(MinIdx+.02,Minima,num2str((1:numel(Minima))'))