% clear all;
% close all;
[openmkfile,openmkpath]=uigetfile('*.smr','Please select the Spike file');
fid=fopen([openmkpath openmkfile],'rb'); %,'ieee-le')
% [WholeLFP,LFPheader]=SONGetChannel(fid, 1);
[WholeLFP,LFPheader]=SONGetChannel(fid, 1);
[WholeEKG,EKGheader]=SONGetChannel(fid, 2);
fclose(fid);
cd(openmkpath)
WholeLFPraw=WholeLFP;
Conv2mV = 5/double(LFPheader.max); % conveting yaxis  values to mV and 5 =max value in Spike2 software
WholeLFP=double(WholeLFP)*Conv2mV;
Fs=1e4;j=(1:length(WholeLFP))/Fs;
figure; plot(j,WholeLFP);

Ictal_st =650; % 1300;
Ictal_ed =2000;%2650;
 
B_st = 770; % 1340;  % second;
B_ed =780;%1360;  % second

LFP = (WholeLFP(Ictal_st*10^4:Ictal_ed*10^4));
% EKG = double(WholeEKG(Ictal_st*10^4:Ictal_ed*10^4));
Fs=1e4;j=(1:length(LFP))/Fs;
figure; plot(j,LFP); 
ylabel ('mV');xlabel('sec');axis([-inf inf -inf inf]); title('raw')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select portion
ephys1=LFP;
% ephys1=ephys(1:8e05);
% [f fftdata]=fftshow(ephys1,10000);suptitle( 'raw')
 
% 1Hz to 50Hz
[c k] = cheby1(2,0.5,[1 50]/5000,'bandpass'); 
q = filtfilt(c,k,ephys1);
[f1 fftdata1]=fftshow(q,10000);
suptitle( 'with 1Hz to 50Hz')
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
figure;
subplot(211);plot(j,ephys1)
subplot(212);plot(j,q2)
LFP=q2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% 
% figure;
% ax1=subplot(211);plot(EKG);title('in optic tectum')
% ax2= subplot(212);plot(LFP);title('in telencephalon')
% linkaxes([ax2,ax1],'x');
% suptitle('ptz  ')

%%%%%%%%%%%%%%%%%



R = 200;
B = ones(R,1);

x = LFP-mean(LFP);
X = x.*x;
X1 = conv(X,R);
LL_LFP = X1(R:length(X1)-R+1);

% detect onset base on line length
baseline = mean(LL_LFP((B_st-Ictal_st)*10000:(B_ed-Ictal_st)*10000)); %baseline LL
b_std = std(LL_LFP((B_st-Ictal_st)*10000:(B_ed-Ictal_st)*10000)); % baseline LL SD
P_Thre = std(LFP((B_st-Ictal_st)*10000:(B_ed-Ictal_st)*10000))*10; % peak threthold
Thre = baseline+b_std; % LL threshold
% ThreFactor=5;

% detect on/off
P_smooth = smooth(LL_LFP,100);
figure;plot(P_smooth)
h=hline(Thre,'m');
h=hline(baseline,'k');

On_OFF = diff(P_smooth>Thre);
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

I = find(Interval < 5000); 
% if the interval between 2 events is less than 0.5 sec, it will be treated as on event

S_start(I+1) = NaN;
S_end(I) = NaN;
S_start(isnan(S_start)) = [];
S_end(isnan(S_end)) = [];

Duration = S_end-S_start;
%% 
for i = 1:length(S_start)
    p1 = LFP(S_start(i):S_end(i));
    p = max(abs(p1));    
    if Duration(i) < 2000 | p < P_Thre 
% if the duration of an event is less than 0.2 s ....
%....and the maximal peak is less than the threthold
        S_start(i) = NaN;
        S_end(i) = NaN;
        Duration(i) = NaN;
    end
    
end

S_start(isnan(S_start)) = [];
S_end(isnan(S_end)) = [];
Duration(isnan(Duration)) = [];
LL_start = S_start;
LL_end = S_end;
LL_Duration = Duration;

%% 
figure
plot(LFP)
hold on
for i = 1:length(LL_start);
    plot(LL_start(i):LL_end(i),LFP(LL_start(i):LL_end(i)),'r')
end
% title('LFP');

figure
subplot(211)
plot(LFP)
hold on
for i = 1:length(LL_start);
    plot(LL_start(i):LL_end(i),LFP(LL_start(i):LL_end(i)),'r')
end
title('LFP');

subplot(212)
plot(LL_LFP);
hold on
for i = 1:length(LL_start);
    plot(LL_start(i):LL_end(i),LL_LFP(LL_start(i):LL_end(i)),'r')
end

title('Line Length')
%% 

LL_Duration = LL_Duration/10000;
figure
bar(LL_Duration)
title('Duration')
hold on
plot([1 length(LL_Duration)],[mean(LL_Duration) mean(LL_Duration)]);


LL_Duration = Duration;
% subplot(121)
hold on;
[cidx,ctrs,SUMD, D] = kmeans(LL_Duration,2);
[Q I] = sort(ctrs);
Ictal = LL_Duration(cidx == I(2));
IIS = LL_Duration(cidx == I(1));

figure; hold on
plot(ones(size(Ictal)),Ictal,'g.','MarkerSize',8)
plot(ones(size(IIS)),IIS,'r.','MarkerSize',8)
plot([1 1],ctrs,'ko','MarkerSize',8);
mean(Q)

% subplot(122)
hold on
[cidx3,ctrs3,SUMD3,D3] = kmeans(LL_Duration,3);
[Q I] = sort(ctrs3);
IIS3 = LL_Duration(cidx3 == I(1));
MED3 = LL_Duration(cidx3 == I(2));
Ictal3 = LL_Duration(cidx3 == I(3));
plot(ones(size(IIS3))*2,IIS3,'r.','MarkerSize',8)
plot(ones(size(MED3))*2,MED3,'b.','MarkerSize',8)
plot(ones(size(Ictal3))*2,Ictal3,'g.','MarkerSize',8)
plot([1 1 1]*2,ctrs3,'ko','MarkerSize',8);
mean(Q(1:2))
mean(Q(2:3))


axis([0 3 0 ceil(max(LL_Duration))])
title('kmean')
ylabel('Duration')

% power in 5 and 50 Hz
Baseline = LFP((B_st-Ictal_st)*10000:(B_ed-Ictal_st)*10000);
[s f t] = spectrogram(Baseline, 1000, 500, 8192,10000);

B_Power = abs(s(1:42,:)); % f(42)=50Hz  f(5)=4.8Hz???
for i = 1:length(S_start)
    p1 = LFP(S_start(i):S_end(i));
    [s f t] = spectrogram(p1, 1000, 500, 8192,10000);
    S_power = abs(s(1:42,:));
   Hz50(i)= sum(S_power(42,:))
   Hz5(i)=sum(S_power(5,:))
    P_A_ratio(i) = sum(S_power(42,:))/sum(S_power(5,:));
end

figure
plot(LL_Duration,1./P_A_ratio,'.');
xlabel('Duration (s)');
ylabel('P50/P5')

Fcount=[mean(Hz50),mean(Hz5), mean(P_A_ratio)];
    figure; bar(Fcount)

%% Number of events Vs Duration histogram
figure;
G = fix(LL_Duration/2000);
N = zeros(1,max(G));
for i = 1:length(LL_Duration)
    N(G(i)) = N(G(i))+1;
end
I = (2:length(N)+1)*0.2;
bar(I,N);
xlabel('Duration (s)');
ylabel('Number')

%% 
figure
plot(LL_Duration,P_A_ratio,'.');
xlabel('Duration (s)');
ylabel('P50/P5')

t = [Ictal_st Ictal_ed B_st B_ed];

% 
% sname = [openmkfile(1:length(openmkfile)-4) '_' num2str(Ictal_st) '_' num2str(Ictal_ed) '.mat'];
% save(sname,'LFP','Duration','P_A_ratio','t');

% plot 3 long seizures

I = find(LL_Duration>4,3);

figure
subplot(2,2,1)
plot(LL_Duration,P_A_ratio,'.');
xlabel('Duration (s)');
ylabel('P50/P5')
hold on;

plot(LL_Duration(I(1:3)),P_A_ratio(I(1:3)),'rs');
hold off;
% I = find(LL_Duration>4,4);
for i = 1:3
    p1 = LFP(S_start(I(i))-10000:S_end(I(i))+10000);
    [s f t] = spectrogram(p1, 1000, 500, 8192,10000);
    S_power = abs(s(1:42,:));
    k = (i>1.5)*2+i+1;
    subplot(4,2,k)
    plot(p1)
    title([num2str(I(i)) '----' num2str(LL_Duration(I(i)))]);
    subplot(4,2,k+2)
    imagesc('XData',t,'YData',f(1:42),'CData',S_power)
    title([num2str(I(i)) '----' num2str(P_A_ratio(I(i)))]);



end

%%  ALL FIGURES IN SUBPLOT 
% close all
f=figure;set(f, 'Position', [4 541 1460 257]);
subplot(141)
Fs=1e4;j=(1:length(LFP))/Fs;
plot(j,LFP); xlabel('sec')
hold on
for i = 1:length(LL_start);
    plot((LL_start(i):LL_end(i))/Fs,LFP(LL_start(i):LL_end(i)),'r')
end
% title('LFP');
axis([-inf inf -inf inf]);

subplot(142)
G = fix(LL_Duration/2000);
N = zeros(1,max(G));
for i = 1:length(LL_Duration)
    N(G(i)) = N(G(i))+1;
end
I = (2:length(N)+1)*0.2;
bar(I,N);
xlabel('Duration (s)');
ylabel('Number')

subplot(143)
hold on
plot(ones(size(Ictal)),Ictal,'g.','MarkerSize',8)
plot(ones(size(IIS)),IIS,'r.','MarkerSize',8)
plot([1 1],ctrs,'ko','MarkerSize',8);
plot(ones(size(IIS3))*2,IIS3,'r.','MarkerSize',8)
plot(ones(size(MED3))*2,MED3,'b.','MarkerSize',8)
plot(ones(size(Ictal3))*2,Ictal3,'g.','MarkerSize',8)
plot([1 1 1]*2,ctrs3,'ko','MarkerSize',8);
axis([0 3 0 ceil(max(LL_Duration))])
title('kmean')
ylabel('Duration')

subplot(144)
plot(LL_Duration,P_A_ratio,'.');
xlabel('Duration (s)');
ylabel('P50/P5')

%%  counts and duration

%%%%% number of events
G = fix(LL_Duration/2000);
N = zeros(1,max(G));
for i = 1:length(LL_Duration)
    N(G(i)) = N(G(i))+1;
end
I = (2:length(N)+1)*0.2;
% bar(I,N);
% xlabel('Duration (s)');
% ylabel('Number')

IIC=N(find(I<=2));
IC=N(find(I>2));
Counts=[sum(IIC), sum(IC), sum(IIC)/sum(IC)]

%%%%% duration
LL_Duration
dIIC=find(LL_Duration/Fs<=2); 
dIC=find(LL_Duration/Fs>2)

if isempty(dIIC); 
    dIICmean =0; 
else
    dIICmean=mean(LL_Duration(dIIC))/Fs
end
dICmean=mean(LL_Duration(dIC))/Fs

%%%%% amplitude
for i = 1:length(LL_start);
%    figure
   spike=LFP(LL_start(i):LL_end(i));
 maxAmp(i)=max(spike)
  end
figure;bar(maxAmp)

%%%%% frequency p50




%%%%%%%  figure
Counts=[sum(IIC), sum(IC), sum(IIC)/sum(IC); dIICmean, dICmean,dIICmean/dICmean ]
figure
if sum(IIC)==0
    bar3(Counts)
else
    bar3(length(Counts),Counts)
end
set(gca,'XTickLabel',{'IIC','IC','IIC/IC'})
set(gca,'yTickLabel',{'Count','Duration'})
labels = arrayfun(@(value) num2str(value,'%2.1f'),Counts,'UniformOutput',false);
for ii = 1:size(d,1)
  for jj = 1:size(d,2)
    TH=text(jj,ii,d(ii,jj),num2str(d(ii,jj)))
    set(TH,'Horizontalalignment','center',...
'verticalalignment','bottom') ;
  end
end
  


  

%%%%%%%%%%%%%%%
figure;
I = find(LL_Duration<2,3);
subplot(2,2,1)
plot(LL_Duration,P_A_ratio,'.');
xlabel('Duration (s)');
ylabel('P50/P5')
hold on;
plot(LL_Duration(I(1:3)),P_A_ratio(I(1:3)),'rs');
hold off;

for i = 1:3
    p1 = LFP(S_start(I(i))-10000:S_end(I(i))+10000);
    [s f t] = spectrogram(p1, 1000, 500, 8192,10000);
    S_power = abs(s(1:42,:));
    P_A = sum(S_power(42,20:length(t)-20))/sum(S_power(5,20:length(t)-20));
    k = (i>1.5)*2+i+1;
    subplot(4,2,k)
    plot(p1)
    title([num2str(I(i)) '----' num2str(LL_Duration(I(i)))]);
    subplot(4,2,k+2)
    imagesc('XData',t,'YData',f(1:42),'CData',S_power)
    title([num2str(I(i)) '----' num2str(P_A_ratio(I(i)))]);
end






