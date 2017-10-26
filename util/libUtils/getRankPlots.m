function f= getRankPlots(ranks,thresh)
% script to plot Ranks on electrode grid, assuming a square grid with electrode one  at
% bottom right and last electrode  at top left
%
%inputs
% ranks:  of electrodes
% thresh: num of top electrodes  to show, default: shows all

if nargin<2
    thresh=[];
end


numElectrodes = length(ranks);
gridSide = sqrt(numElectrodes);

matRank = rot90(reshape((ranks<thresh).*ranks,gridSide,gridSide),2);
channelMatrix = matRank;
% channelMatrix = rot90(reshape(1:numElectrodes,gridSide,gridSide),2); % to
% show channel numbers
textstr = num2str(channelMatrix(:),'%d');
textstr = strtrim(cellstr(textstr));

% set color of text str here
textColors = repmat([1,0,0],256,1);

f=imagesc(matPrank);colormap gray
[x,y] = meshgrid(1:gridSide);
hStr = text(x(:),y(:),textstr(:),'HorizontalAlignment','center');
set(hStr,{'Color'},num2cell(textColors,2));
% set(gca,'xtick',[],'ytick',[])
%
set(gca,'xaxisLocation','top','yaxisLocation','right')

set(gca,'xtick',1:gridSide,'XTickLabel',[numElectrodes:-gridSide:gridSide])
set(gca,'ytick',1:gridSide,'YTickLabel',[gridSide:-1:1])

end