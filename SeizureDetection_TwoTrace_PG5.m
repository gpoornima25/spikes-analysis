clear all; 
close all;
[openmkfile,openmkpath]=uigetfile('*.smr','Please select the Spike file');
fid=fopen([openmkpath openmkfile],'rb'); %,'ieee-le')
[WholeLFP,LFPheader]=SONGetChannel(fid, 1);
[WholeEKG,EKGheader]=SONGetChannel(fid, 2);
fclose(fid);cd(openmkpath);

Ictal_st = 1800;
Ictal_ed = 4200;
CutOff=0.7;

LFP = double(WholeLFP(Ictal_st*10^4:Ictal_ed*10^4));
EKG = double(WholeEKG(Ictal_st*10^4:Ictal_ed*10^4));%  EKG=EKG/10;% gain was off for first few tel-ot exps
% LFP=LFP*1.2;
%   EKG=EKG/20;
%  EKG=EKG*1.3;

B_LFP = LFP; B_EKG = EKG;
R = 200;
B = ones(R,1);
[c f] = butter(1,[1/5000 50/5000], 'bandpass');
LFP = filtfilt(c,f,LFP);
X = LFP.*LFP;
X1 = conv(X,R);

EKG = filtfilt(c,f,EKG);
LL_LFP = X1(R:length(X1)-R+1);
X = EKG.*EKG;
X1 = conv(X,R);
LL_EKG = X1(R:length(X1)-R+1);

figure;
Fs=1e4;j=(1:length(LFP))/Fs;
subplot(211);plot(j,LFP);axis([-inf inf -inf inf])
subplot(212);plot(j,EKG);axis([-inf inf -inf inf])
%%

figure;
ax1=subplot(211);
plot(LL_LFP)
Thre_LFP =1e6;
h=hline(Thre_LFP,'k');


ax2=subplot(212);
plot(LL_EKG)
Thre_EKG =1e8;
h=hline(Thre_EKG,'k');

linkaxes([ax1,ax2],'x')

%% LFP 

% detect on/off of LFP
P_smooth = smooth(LL_LFP,100); LL_smooth_LFP = P_smooth;
On_OFF = diff(P_smooth>Thre_LFP);
S_end = find(On_OFF == -1);
S_start = find(On_OFF == 1);

if S_start(1)>S_end(1)
    S_start = [1 S_start'];
    S_start = S_start';
end

if length(S_start) > length(S_end)
    S_end = [S_end' length(LL_LFP)];
    S_end = S_end';
end

l = length(S_start);

clear Interval;
for i = 1:l-1
    Interval(i) = S_start(i+1) - S_end(i);
end

I = find(Interval < 20000); % if the interval between 2 events is less than 0.5 sec, it will be treated as on event

S_start(I+1) = NaN;
S_end(I) = NaN;
S_start(isnan(S_start)) = [];
S_end(isnan(S_end)) = [];

Duration = S_end-S_start;
% % 
% for i = 1:length(S_start)
%     p1 = LFP(S_start(i):S_end(i));
%     p = max(abs(p1));    
%     % if Duration(i) < 2000 | p < P_Thre % if the duration of an event is less than 0.2 s and the maximal peak is less than the threthold
%     if Duration(i) < 2000
%         S_start(i) = NaN;
%         S_end(i) = NaN;
%         Duration(i) = NaN;
%     end
%     
% end

S_start(isnan(S_start)) = [];
S_end(isnan(S_end)) = [];
Duration(isnan(Duration)) = [];
LFP_start = S_start;
LFP_end = S_end;
LFP_Duration = Duration;
% 
% f3=figure;
% subplot(211)
% plot(B_LFP)
% hold on
% for i = 1:length(LFP_start);
%     plot(LFP_start(i):LFP_end(i),B_LFP(LFP_start(i):LFP_end(i)),'r')
% end
% title('LFP');
% 
% subplot(212)
% plot(LL_LFP);
% hold on
% for i = 1:length(LFP_start);
%     plot(LFP_start(i):LFP_end(i),LL_LFP(LFP_start(i):LFP_end(i)),'r')
% end
% title('Line Length')
% saveppt('051217.ppt','TEL_LineLength','-f3')


% Cross corr between LFP and EKG spikes
[c f] = butter(1,1.5/5000,'high');

LFP = filtfilt(c,f,LFP);
EKG = filtfilt(c,f,EKG);
j=1;Fs=10^4;
for i = 1:length(LFP_start)-1
    L1 = LFP(LFP_start(i)-5000:LFP_end(i)+5000);
    L2 = EKG(LFP_start(i)-5000:LFP_end(i)+5000);
    [r,p] = corrcoef(L1,L2);
    R(j) = r(1,2);
    [C,Lags]=xcorr(L1,L2);
    [~,M] = max(abs(C));
    SampleDiff = Lags(M); % t21 = finddelay(T1,S)
    timeDiff = round((SampleDiff/Fs)*1e3);
    T(j)=timeDiff;
    j=j+1;
end

I1 = find(R<CutOff);
I=I1;
f4=figure;
set(f4, 'Position', [1 41 1600 783]);
subplot(3,4,1)
bar(R)
hold on
plot(I,R(I),'r*')
axis([-inf inf -inf inf]);
j=1;
for i = 1:min([11 length(I)])
    L1 = LFP(LFP_start(I(i))-5000:LFP_end(I(i))+5000);
    L2 = EKG(LFP_start(I(i))-5000:LFP_end(I(i))+5000);
    subplot(3,4,j+1);j=j+1;
    plot(L1); axis([-inf inf -inf inf])
    hold on;
    plot(L2);axis([-inf inf -inf inf])
    title([num2str(I(i)) ' (' num2str(round(R(I(i)),2))  ') OT leads hb by ' num2str(T(I(i))) 's'] ,'FontSize', 6)
end
legend('OT', 'Hb')
tightfig
suptitle('OT')
% saveppt('051217.ppt','OT','-f4')
if i==11
    f=figure;j=1;set(f, 'Position', [1 41 1600 783]);
    for i = 12:length(I)
    L1 = LFP(LFP_start(I(i))- 5000:LFP_end(I(i))+5000);
        L2 = EKG(LFP_start(I(i))-5000:LFP_end(I(i))+5000);
        subplot(5,4,j);j=j+1;
        plot(L1); 
        hold on;
        plot(L2);
        title([num2str(I(i)) ' (' num2str(round(R(I(i)),2))  ') OT leads hb by ' num2str(T(I(i))) 's'],'FontSize', 6)
    end
    tightfig

else    
end

%% EKG
P_smooth = smooth(LL_EKG,100);  LL_smooth_EKG = P_smooth;
On_OFF = diff(P_smooth>Thre_EKG);
S_end = find(On_OFF == -1);
S_start = find(On_OFF == 1);

if S_start(1)>S_end(1)
    S_start = [1 S_start'];
    S_start = S_start';
end

if length(S_start) > length(S_end)
    S_end = [S_end' length(LL_LFP)];
    S_end = S_end';
end

l = length(S_start);

clear Interval;
for i = 1:l-1
    Interval(i) = S_start(i+1) - S_end(i);
end

I = find(Interval < 20000); % if the interval between 2 events is less than 0.5 sec, it will be treated as on event

S_start(I+1) = NaN;
S_end(I) = NaN;
S_start(isnan(S_start)) = [];
S_end(isnan(S_end)) = [];

Duration = S_end-S_start;
% % 
% for i = 1:length(S_start)
%     p1 = LFP(S_start(i):S_end(i));
%     p = max(abs(p1));    
%     % if Duration(i) < 2000 | p < P_Thre % if the duration of an event is less than 0.2 s and the maximal peak is less than the threthold
%     if Duration(i) < 2000
%         S_start(i) = NaN;
%         S_end(i) = NaN;
%         Duration(i) = NaN;
%     end
%     
% end

S_start(isnan(S_start)) = [];
S_end(isnan(S_end)) = [];
Duration(isnan(Duration)) = [];
EKG_start = S_start;
EKG_end = S_end;
EKG_Duration = Duration;

%%%%%%
f3=figure;
Fs=1e4;j=(1:length(B_LFP))/Fs;
subplot(211);plot(j,B_LFP);xlabel('sec')

hold on
for i = 1:length(LFP_start)
    k=(LFP_start(i):LFP_end(i))/Fs;
    plot(k,B_LFP(LFP_start(i):LFP_end(i)),'r')
%     plot((1:length(B_LFP(LFP_start(i):LFP_end(i))))/Fs,B_LFP(LFP_start(i):LFP_end(i)),'r')

end
title('OT'); axis([-inf inf min(B_LFP) max(B_LFP) ]);

Fs=1e4;j=(1:length(B_EKG))/Fs;
subplot(212);plot(j,B_EKG);xlabel('sec')
hold on
for i = 1:length(EKG_start);
        k=(EKG_start(i):EKG_end(i))/Fs;
    plot(k,B_EKG(EKG_start(i):EKG_end(i)),'r')
end
title('HB'); axis([-inf inf min(B_LFP) max(B_LFP) ]);

% saveppt('051217.ppt','OT_LineLength','-f5')
%%%%%%%%%%%

[c f] = butter(1,1.5/5000,'high');
LFP = filtfilt(c,f,LFP);
EKG = filtfilt(c,f,EKG);
j=1;Fs=10^4;
for i = 2:length(EKG_start)-1
    L1 = LFP(EKG_start(i)-5000:EKG_end(i)+5000);
    L2 = EKG(EKG_start(i)-5000:EKG_end(i)+5000);
    [r,p] = corrcoef(L1,L2);
    R_EKG(j) = r(1,2);
    [C,Lags]=xcorr(L1,L2);
    [~,M] = max(abs(C));
    SampleDiff = Lags(M); % t21 = finddelay(T1,S)
    timeDiff = round((SampleDiff/Fs)*1e3);
    T_EKG(j)=timeDiff;
     j=j+1;
end

I2 = find(R_EKG<CutOff);
clear I
I=I2;
f6=figure;
set(f6, 'Position', [1 41 1600 783]);
subplot(3,4,1)
bar(R_EKG)
hold on
plot(I,R_EKG(I),'r*')
axis([-inf inf -inf inf]);
j=1;
for i = 1:min([11 length(I)])
    L1 = LFP(EKG_start(I(i))-5000:EKG_end(I(i))+5000);
    L2 = EKG(EKG_start(I(i))-5000:EKG_end(I(i))+5000);
    subplot(3,4,j+1);j=j+1;
    plot(L1); axis([-inf inf -inf inf])
    hold on;
    plot(L2);axis([-inf inf -inf inf])
    title([num2str(I(i)) ' ('  num2str(round(R_EKG(I(i)),2)) ') OT leads hb by ' num2str(T_EKG(I(i))) 's'] ,'FontSize', 6)
end
legend('OT', 'HB')
tightfig;
suptitle('HB')
% saveppt('051217.ppt','OT spikes against OT','-f6')

if i==11
    f7=figure;
    set(f7, 'Position', [1 41 1600 783]);j=1;
    for i = 12:length(I)% adjust this number based on length of I
        L1 = LFP(EKG_start(I(i))-5000:EKG_end(I(i))+5000);
        L2 = EKG(EKG_start(I(i))-5000:EKG_end(I(i))+5000);
        subplot(5,4,j)
        plot(L1); hold on;
        plot(L2);j=j+1;
    title([num2str(I(i)) ' (' num2str(round(R_EKG(I(i)),2)) ') OT leads hb by ' num2str(T_EKG(I(i))) 's'] ,'FontSize', 6)
    end
else
end

clc
size(R)
size(R_EKG)
size(I1)
size(I2)
mean(R)
mean(R_EKG)

% saveppt('051217.ppt','contd...','-f')

% 
% I = find(R<CutOff);
% f=figure;
% set(f, 'Position', [1 41 1600 783]);
% Fs=10^4;
% for i = 1:length(I)
%     L1 = LFP(LFP_start(I(i))-25000:LFP_end(I(i))+25000);
%     L2 = EKG(LFP_start(I(i))-25000:LFP_end(I(i))+25000);
% [C2,Lags2]=xcorr(L1,L2);
% [C2lfp,Lags2lfp]=xcorr(L1,L1);
% [C2ekg,Lags2ekg]=xcorr(L2,L2);
%     subplot(5,5,i)
% h1=plot(Lags2/Fs,C2);hold on
% h2=plot(Lags2lfp/Fs,C2lfp,'r'); 
% h3=plot(Lags2ekg/Fs,C2ekg,'g'); 
% h=vline(0,'y');
% ylabel('correlation');xlabel('sec');
% axis([-1 1 -inf inf])
% [~,M] = max(abs(C2));
% SampleDiff = Lags2(M); % t21 = finddelay(T1,S)
% timeDiff = round((SampleDiff/Fs)*1e3);
% title([num2str(I(i)) ') tel leads ot by ' num2str(timeDiff) 's'],'FontSize', 6)
% end
% % legend([h1 h2 h3],{'tel ot', 'tel', 'ot'})
% tightfig
% saveppt('051217.ppt','crosscorr of TEL spikes','-f')

% saveppt('051217.ppt','contd...','-f7')

% 
% f8=figure;
% set(f8, 'Position', [1 41 1600 783]);
% Fs=10^4;
% for i = 1:length(I)
%     L1 = LFP(EKG_start(I(i))-25000:EKG_end(I(i))+25000);
%     L2 = EKG(EKG_start(I(i))-25000:EKG_end(I(i))+25000);
% [C2,Lags2]=xcorr(L1,L2);
% [C2lfp,Lags2lfp]=xcorr(L1,L1);
% [C2ekg,Lags2ekg]=xcorr(L2,L2);
%     subplot(6,6,i)
% h1=plot(Lags2/Fs,C2);hold on
% h2=plot(Lags2lfp/Fs,C2lfp,'r'); 
% h3=plot(Lags2ekg/Fs,C2ekg,'g'); 
% h=vline(0,'y');
% ylabel('correlation');xlabel('sec');
% axis([-1 1 -inf inf])
% [~,M] = max(abs(C2));
% SampleDiff = Lags2(M); % t21 = finddelay(T1,S)
% timeDiff = round((SampleDiff/Fs)*1e3);
% title([num2str(I(i)) ') tel leads ot by ' num2str(timeDiff) 's'],'FontSize', 6)
% end
% legend([h1 h2 h3],{'tel ot', 'tel', 'ot'})
% tightfig;
% % saveppt('051217.ppt','crosscorr of OT spikes','-f8')
