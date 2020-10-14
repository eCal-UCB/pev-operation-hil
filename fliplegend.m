function []=fliplegend( labels )
fig=gcf;
bb=fig.Children.Children';
l=legend(bb, labels,'Location','Best');
legholder=l.String;
for i=1:length(legholder)
    l.String{i}=legholder{end+1-i};
end
end