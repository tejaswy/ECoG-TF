function plotChannelgrid(Channels,titleStr)


channelMatrix = rot90(reshape(1:256,16,16),2);
textstr = num2str(channelMatrix(:),'%d');
textstr = strtrim(cellstr(textstr));
textColors = repmat([1,0,0],256,1);

p=zeros(256,1);
p(Channels)=1;
matPrank = rot90(reshape(p,16,16),2);
imagesc(matPrank);colormap gray; colorbar
[x,y] = meshgrid(1:16);
hStr = text(x(:),y(:),textstr(:),'HorizontalAlignment','center');
set(hStr,{'Color'},num2cell(textColors,2));
set(gca,'xtick',[],'ytick',[])
title(titleStr)

end
