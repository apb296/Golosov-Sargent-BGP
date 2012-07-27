function [TGrid,X]=getPlotData(dd)
TGridSize=10*length(dd);
TGrid=[];
X=[];
for i=1:length(dd)
TGrid=[TGrid linspace(i-1,i,TGridSize/length(dd))];
X=[X ones(1,TGridSize/length(dd))*dd(i)];
end

end