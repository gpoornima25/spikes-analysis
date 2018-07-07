clear all;
 close all;
[openmkfile,openmkpath]=uigetfile('*.smr','Please select the Spike file');
fid=fopen([openmkpath openmkfile],'rb'); %,'ieee-le')
% [WholeLFP,LFPheader]=SONGetChannel(fid, 1);
[WholeLFP,LFPheader]=SONGetChannel(fid, 1);
[WholeEKG,EKGheader]=SONGetChannel(fid, 2);
% figure; plot(WholeEKG)
fclose(fid);
cd(openmkpath)



Ictal_st = 280;
Ictal_ed = 2460;
 
B_st = 2360;  % second;
B_ed = 2450;  % second




Fs=10^4;

LFP = double(WholeLFP(Ictal_st*10^4:Ictal_ed*10^4));
EKG = double(WholeEKG(Ictal_st*10^4:Ictal_ed*10^4));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select portion
ephys1=LFP;
% ephys1=ephys(1:8e05);
% [f fftdata]=fftshow(ephys1,10000);suptitle( 'raw')
 
% 1Hz to 55Hz
[c k] = cheby1(2,0.5,[1 55]/5000,'bandpass'); 
q = filtfilt(c,k,ephys1);
[f1 fftdata1]=fftshow(q,10000);suptitle( 'with 1Hz to 55Hz')
close

%notch filter specs
Fs=10000; t = (0:length(q)-1)/Fs;
%2Hz notch
    d = designfilt('bandstopiir','FilterOrder',2, ...
                   'HalfPowerFrequency1',1,'HalfPowerFrequency2',3, ...
                   'DesignMethod','butter','SampleRate',Fs);
%     fvtool(d,'Fs',Fs)
    q2 = filtfilt(d,q);
[f2 fftdata2]=fftshow(q2,10000);suptitle( ' 2Hz notch filter ')
close
% figure;
% subplot(211);plot(ephys1)
% subplot(212);plot(q2)
LFP=q2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select portion
ephys1=EKG;
% ephys1=ephys(1:8e05);
% [f fftdata]=fftshow(ephys1,10000);suptitle( 'raw')
 
% 1Hz to 55Hz
[c k] = cheby1(2,0.5,[1 55]/5000,'bandpass'); 
q = filtfilt(c,k,ephys1);
[f1 fftdata1]=fftshow(q,10000);suptitle( 'with 1Hz to 55Hz')
close

%notch filter specs
Fs=10000; t = (0:length(q)-1)/Fs;
%2Hz notch
    d = designfilt('bandstopiir','FilterOrder',2, ...
                   'HalfPowerFrequency1',1,'HalfPowerFrequency2',3, ...
                   'DesignMethod','butter','SampleRate',Fs);
%     fvtool(d,'Fs',Fs)
    q2 = filtfilt(d,q);
[f2 fftdata2]=fftshow(q2,10000);suptitle( ' 2Hz notch filter ')
close
% figure;
% subplot(211);plot(ephys1)
% subplot(212);plot(q2)
EKG=q2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% figure;
% ax1=subplot(211);plot(EKG-mean(EKG));title('in optic tectum');axis([-inf inf -inf inf])
% ax2= subplot(212);plot(LFP-mean(LFP));title('in telencephalon');axis([-inf inf -inf inf])
% linkaxes([ax2,ax1],'x');
% suptitle('ptz  ')

%% in Secs
TotTime=length(LFP)/1e4;


%% cross correlation
 [C2,Lags2]=xcorr(LFP,EKG);
[C2lfp,Lags2lfp]=xcorr(LFP,LFP);
[C2ekg,Lags2ekg]=xcorr(EKG,EKG);

% hFig=figure;subplot(121)
% plot(Lags2/.1e4,C2*12); hold on;
% plot(Lags2lfp/.1e4,C2lfp*100,'r'); 
% plot(Lags2ekg/.1e4,C2ekg/2,'g');
% ylabel('correlation')
% xlim([-10 10]);xlabel('sec') % 3.2e7=3200sec of recording
% legend('tel hb', 'tel', 'hb')


figure
h1=plot(Lags2/Fs,C2);hold on
h2=plot(Lags2lfp/Fs,C2lfp,'r'); 
h3=plot(Lags2ekg/Fs,C2ekg,'g'); 
h=vline(0,'y');
ylabel('correlation')
xlim([-.1 .1]);xlabel('sec') %
suptitle('ot-hb #5')
% line([0 0], [0 inf],'Color',[0.5 0.5 0.5],'LineStyle','--')
% ylim([3e10 12e10])
legend([h1 h2 h3],{'ot hb', 'ot', 'hb'})
ylim([0 4e12])
 print('-clipboard','-dmeta')
% hgexport(hFig,'-clipboard')

% selcting subplot
% hd = findall(gcf,'type', 'axes'); %

%% detecting time diff

% T2=LFP;%template
% S=EKG;%signal

[C2,Lags2] = xcorr(LFP,EKG);        
Fs=10^4;
% figure
% ax(1) = subplot(2,1,1); 
% plot(lag1/Fs,C1,'k')
% ylabel('Amplitude')
% grid on
% title('Cross-correlation between T1 and S')

% high peak in the second subplot indicates that
% signal is present in the  template.

[~,I] = max(abs(C2));
SampleDiff = Lags2(I) % t21 = finddelay(T1,S)
timeDiff = SampleDiff/Fs % round(Timediff *1e3)

% The peak of the cross correlation implies
% that the signal is present in template T2 (tel)
% starting after 'round(Timediff *1e3)' ms. 
% Ex: tel-ot#2LFP->tel leads EKG->ot by 2ms 
% tel-ot#3LFP->tel leads EKG->ot by .2ms 

% In other words, 
% signal T1 leads signal S by '(SampleDiff)' samples as
% indicated by SampleDiff. This information 
% can be used to align the signals



