function  videoMatrix= visualizeActivity(data,Fs,toNorm,fps, movingwin , toSave, behaviorSeries)
% Visualize activity in a video
% Inputs:
% --required:
% data : Activity to plot, currently assumes data is from a square grid
%   columns of data are assumed to be time samples. rows could be channels
%   or different trials
%  Fs:  sampling frequency of data
% --optional
% toNorm:=1 to zscore the data along rows. So each column is normalized
%     Default=1
%  fps: frames per second in the video: default: 1
%  movingwin=  (in the form [window winstep] i.e length of moving
%                                              window and step size in seconds)
%           Default: [0.1,0.025] i.e. 100 ms and 25 ms
%  So by default you will see 100 ms averaged activity in 1 frame per
%  second
% toSave  = 1, saves the video in the same folder with date, time as file name
%         =0 ; if you dont want to save. Default=0;
%
%behaviorSeries = num behavioral Sensors x timesamples.. time series with
%behavioral data with same Fs as neural data
% Outputs:
%  videoMatrix


if nargin<2
    error('data and Fs are required')
end

if nargin<3||isempty(toNorm)
    toNorm =1;
end

if nargin<4||isempty(fps)
    fps=1;
end

if nargin<5||isempty(movingwin)
    movingwin=[0.1,0.025];
end

if nargin<6||isempty(toSave)
    toSave=0;
end

if nargin<7|| isempty(behaviorSeries)
    isSubplot=0;
else
    isSubplot=1;
end

nRows= size(data,1);

H= sqrt(nRows);% currently assuming data is from square grid
W=H;
Nwin=round(Fs*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*Fs); % number of samples to step through
N= size(data,2);
winstart=1:Nstep:N-Nwin+1;
nFrames=length(winstart);

if toNorm
    data= zscore(data,0,2);
end


videoMatrix = zeros(H,W,nFrames);

% First put all frames in a 3D matrix
for i=1:nFrames
    temp = mean(data(:,winstart(i):winstart(i)+Nwin-1),2);
    videoMatrix(:,:,i)= reshape(temp,H,H);
    clear temp
end

% text in the boxes has row number in data matrix
channelMatrix = reshape(1:nRows,H,H);
textstr = num2str(channelMatrix(:),'%d');
textstr = strtrim(cellstr(textstr));
textColors = repmat([1,0,0],nRows,1);
[x,y] = meshgrid(1:H);

hStr = text(x(:),y(:),textstr(:),'HorizontalAlignment','center');
f = figure;  colormap gray;
axis tight manual
ax = gca;
% ax.NextPlot = 'replaceChildren';
if isSubplot
    
    for i=1:nFrames
        subplot(211)
        imagesc(squeeze(videoMatrix(:,:,i)));   %colorbar
        hStr = text(x(:),y(:),textstr(:),'HorizontalAlignment','center');
        set(hStr,{'Color'},num2cell(textColors,2));
        set(gca,'xtick',[],'ytick',[])
        tstart= (winstart(i)/Fs)*1000;
        tstop= ((winstart(i)+Nwin-1)/Fs)*1000;
        title(strcat('time in milli sec:' ,num2str(tstart),' to ',num2str(tstop)));
        
        subplot(212)
        plot(behaviorSeries(:,winstart(i):winstart(i)+Nwin-1))
        pause(1/fps);
        
    end
    imagesc(squeeze(videoMatrix(:,:,i)));   %colorbar
    hStr = text(x(:),y(:),textstr(:),'HorizontalAlignment','center');
    set(hStr,{'Color'},num2cell(textColors,2));
    set(gca,'xtick',[],'ytick',[])
    tstart= (winstart(i)/Fs)*1000;
    tstop= ((winstart(i)+Nwin-1)/Fs)*1000;
    title(strcat('time in milli sec:' ,num2str(tstart),' to ',num2str(tstop)));
    
else
    
    
end
if toSave
    
    %Build a name for the matfile
    datetime=datestr(now);
    datetime=strrep(datetime,':','_'); %Replace colon with underscore
    datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
    datetime=strrep(datetime,' ','_');%Replace space with underscore
    filename= strcat(datetime,'videoMatrix');
    fprintf('\n Saving videoMatrix in the current folder')
    save(fileName,'videoMatrix');
end
end



