clc
% clear all
close all
set(0,'defaulttextinterpreter','tex')
fs=400;
movingwin         = [0.25 0.05]; % moving window size in [s]
params.Fs         = fs; % sampling frequency
params.fpass      = [0 200]; % frequencies of interest
params.tapers     = [4 7];
params.trialave   = 1;
p                 = 0.05;
params.err        = [2 p] ;
dname = '/home/tejaswy/ecog_analysis/ecog_analysis/ecog_analysis';
dnamefeat = '/home/tejaswy/ecog_analysis/ecog_analysis/ecog_analysis/NologFeatures';
addpath(genpath('/home/tejaswy/ecog_analysis/ecog_analysis/chronux_2_10'))
addpath(genpath('/home/tejaswy/ecog_analysis/ecog_analysis/ecog_analysis/Misc'))
h=figure;
for subject =1
%     h=subplot(1,3,subject);
    %     [sensory, motor,stg,larynx] = regions(subject);
    %     %%%%%%%%%%%%%%Load ecog data and pitch data
    
    if subject == 1
        ecogAddr                  = fullfile(dname,'/CleanData/ec56.mat');
        load(ecogAddr)
        labelfile                 = fullfile(dname,'Labels/ECOG_ind_labels/EC56Labels_ecogInd.mat');
        load(labelfile)
        PNLabel = double(ec56_PNLabel);
        hlmLabel = double(ec56_HLM_Phon_Labels);
        subID = 'EC56';
%         load (fullfile(dname,'SubjectTrials/PN/ec56'));
        load(fullfile(dnamefeat,'BandPowers/ec56_powers'));
    elseif subject == 2
        ecogAddr                  = fullfile(dname,'/CleanData/ec61.mat');
        load(ecogAddr)
        labelfile                 = fullfile(dname,'Labels/ECOG_ind_labels/EC61Labels_ecogInd.mat');
        load(labelfile)
        PNLabel = double(ec61_PNLabel);
        hlmLabel = double(ec61_HLM_Phon_Labels);
        subID = 'EC61';
%         load (fullfile(dname,'SubjectTrials/PN/ec61'));
        load(fullfile(dname,'HilbertPowers/ec61_powers'));
    elseif subject == 3
        ecogAddr                  = fullfile(dname,'/CleanData/ec69.mat');
        load(ecogAddr)
        labelfile                 = fullfile(dname,'Labels/ECOG_ind_labels/EC69Labels_ecogInd.mat');
        load(labelfile)
        PNLabel = double(ec69_PNLabel);
        hlmLabel = double(ec69_HLM_Phon_Labels);
        subID = 'EC69';
%         load (fullfile(dname,'SubjectTrials/PN/ec69'));
        load(fullfile(dname,'HilbertPowers/ec69_powers'));
    end
    [sensory,motor_premotor,STG,larynx,ch]= regions(subject);
%     % el=unique([sensory,motor_premotor,larynx]);
%     el=motor_premotor;
%     % el=1:256;
%     
%     
%     
%     ch=setdiff(el,badchans);
%     r= find(PNLabel(:,3)>200);
%     PNLabel = PNLabel(r,:);
%     % parameters for spectrogram
%     ch=chan(subject,:);
    
%     ch = chAnova.BP(1:5);
%     for i=1:length(ch)
%         if ch(i)>256
%             ch(i)= ch(i)-256;
%         end
%     end
    
    sp = find(PNLabel(:,4)==1);
    
    sil = find(PNLabel(:,4)==0);
    
    spMarkers = [PNLabel(sp,1),PNLabel(sp,2)];
    edata = data(ch,:);
    [ spS, f, Serrsp ]= mtspectrumc_unequal_length_trials( edata',movingwin, params, spMarkers );
    Serrsp = mean(Serrsp,3);
    plot_vector(mean(spS'),f,'l',[],'k'); hold on
    plot(f,10*log10(Serrsp(1,:)),'k:');
    plot(f,10*log10(Serrsp(2,:)),'k:');
    
    silMarkers = [PNLabel(sil,1),PNLabel(sil,2)];
    [ silS, f, Serr ]= mtspectrumc_unequal_length_trials( edata',movingwin, params, silMarkers );
    plot_vector(mean(silS'),f,'l',[],'r');hold on
    Serr = mean(Serr,3);
    
    plot(f,10*log10(Serr(1,:)),'m:');
    plot(f,10*log10(Serr(2,:)),'m:');
    title(subID);
xlabel('Frequency [in Hz]')

ylabel('10*log_{10}(S)','interpreter','tex')
end
ax= get(h,'Children');
legend(ax([1,2]),'Voc','Sil')


% suptitle('Average Power Spectra of selected electrodes')
% plot(1,1,'k');
% plot(0,0,'r')
% legend('Voc','Sil')
%
%
% % power modulations
% % select electrodes top 5 in HP
%
% spInd = [PNLabel(sp,1)-200,PNLabel(sp,1)+200];
% silInd = [PNLabel(sil,1)-200,PNLabel(sil,1)+200];
%
% ch = chAnova.BP(1:5);
%
% for i=1:5
%     if ch(i)>256
%         ch(i)=ch(i)-256;
%         for k=1:length(sp)
%             spAct(k,:)=(HGP(ch(i),spInd(k,1): spInd(k,2)));
%         end
%         for k=1:length(sil)
%             silAct(k,:)=(HGP(ch(i),silInd(k,1): silInd(k,2)));
%         end
%     else
%         for k=1:length(sp)
%             spAct(k,:)=(BP(ch(i),startInd(k,1): startInd(k,2)));
%         end
%         for k=1:length(sil)
%             silAct(k,:)=(BP(ch(i),silInd(k,1): silInd(k,2)));
%         end
%     end
%     meanactSP = tsmovavg(mean(spAct),'s',10);
%     sterSP  = tsmovavg(std(spAct)/sqrt(length(sp)),'s',10);
%     shadedErrorBar(linspace(-500,500,length(meanactSP)),meanactSP,sterSP  ,'r',1);  hold on
%
%     meanactSil = tsmovavg(mean(silAct),'s',10);
%     sterSil  = tsmovavg(std(silAct)/sqrt(length(sil)),'s',10);
%     shadedErrorBar(linspace(-500,500,length(meanactSil)),meanactSil,sterSil  ,'k',1);  hold off
%     clear spAct
%     clear silAct
%     pause
%
% end
% %
% %
% %
% % for ch = larynx
% %     count=1;
% %     for i=1:length(sp)
% %         if length(PNLabel(sp(i),1):PNLabel(sp(i),2))>256
% %         edata =  data(ch,PNLabel(sp(i),1):PNLabel(sp(i),2));
% %         size(edata)
% %         [spS(count,:),f] = mtspectrumc(edata',params);
% %         count=count+1;
% %         size(spS)
% %         end
% %     end
% %     plot_vector(mean(spS),f,'l',[],'r');
% %
% %     for i=1:length(sil)
% %         edata =  data(ch,PNLabel(sil(i),1): PNLabel(sil(i),2));
% %         [silS(i,:),f] = mtspectrumc(edata',params);
% %         size(silS)
% %     end
% %     plot_vector(mean(silS),f,'l',[],'k');
% %
% %     legend('Voc','Sil')
% %     pause
% %     clc
% % end
%
% %
% % ch =larynx;
% % % speaking
% % for i=1:length(sp)
% %     edata =  data(ch,PNLabel(sp(i),1): PNLabel(sp(i),2));
% %     [spS(i,:),f,Serr] = mtspectrumc(edata',params);
% % end
% % plot_vector(mean(spS),f,'l',mean(Serr),'r',2);
% % hold on;
% %
% %
% % for i=1:length(sil)
% %     edata =  data(ch,PNLabel(sil(i),1): PNLabel(sil(i),2));
% %     [nspS(i,:),f,Serr] = mtspectrumc(edata',params);
% % end
% % plot_vector(mean(nspS),f,'l',mean(Serr),'k',2);
% %
% %
