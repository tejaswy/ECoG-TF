function plotCM(CM, toNorm,groupNames)
% plot confusion matrix (grey scale with values in red)
%inputs:
%CM: Confusion matrix
% toNorm=1, normalize CM,
%        =0 otherwise
% groupNames=  cell with string containing group names Eg: {'Voc' , 'Silent'}
% output
% f=figure handle
if nargin<2
    toNorm=0;
end
if nargin<3
    groupNames={};
end

% f=figure;
numCat = length(CM);
if toNorm==1   % Normalize
    CM = CM./repmat(sum(CM,2),1,numCat);
end

% f=figure
imagesc(CM)
colormap(flipud(gray));
textstr = num2str(CM(:),'%0.2f');
textstr = strtrim(cellstr(textstr));
[x,y] = meshgrid(1:numCat);
hStr = text(x(:),y(:),textstr(:),'HorizontalAlignment','center');

textColors = repmat([1,0,0],numCat^2,1);
set(hStr,{'Color'},num2cell(textColors,2));


xlabel('Predicted Label'); %,'FontSize',14,'Fontweight','bold')
ylabel('True Label') %,'FontSize',14,'Fontweight','bold')


if isempty(groupNames)
    for i=1:numCat     
        groupNames(i)={num2str(i)};
    end
end
set(gca,'XTick',1:numCat,...
    'XTickLabel',groupNames,...
    'YTick',1:numCat,...
    'YTickLabel',groupNames)

end