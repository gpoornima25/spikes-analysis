clear all
close all
clc

%% split file into baseline, PTZ nad paralyzing agent (PA) ( PB- pancuronium bromide , AB- alpha bungratoxin)

%from  CutSpikeTrace.m
[openmkfile,openmkpath]=uigetfile('*.smr','Please select the Spike file');
fid=fopen([openmkpath openmkfile],'rb'); %,'ieee-le')
cd(openmkpath)
[WholeLFP,LFPheader]=SONGetChannel(fid, 1);
 [WholeEKG,EKGheader]=SONGetChannel(fid, 2);
 figure; plot(WholeEKG)
% figure;plot(WholeLFP)

% % for file with part2
% [openmkfile,openmkpath]=uigetfile('*.smr','Please select the Spike file');
% fid=fopen([openmkpath openmkfile],'rb'); %,'ieee-le')
% % [WholeLFP,LFPheader]=SONGetChannel(fid, 1);
% [WholeEKG2,EKGheader2]=SONGetChannel(fid, 2);
% figure; plot(WholeEKG2)

%manually selct baseline, post PTZ and post PA(paralyzing agent) parts referring from spike 2
% baseline=WholeEKG(1:9.797e06);
% postPTZ=WholeEKG(1.34e07:1.375e07);
postPA=WholeEKG(3.9e06:length(WholeEKG));
figure; 
subplot(311)
plot (baseline); axis([-inf inf -10000 10000]);title( 'baseline')
subplot(312)
plot(postPTZ); axis([-inf inf  -10000 10000]);title( 'post PTZ')
subplot(313)
plot(postPA); axis([-inf inf -10000 10000]);title( 'post PA')
suptitle('09152016 with pa Fish3')


%% Filtering (noise about 0.05 to 0.1mV)

%select small part from the baseline with no spikes
BaselineSample=baseline(6.9e06:7.4e06);
% BaselineMean= mean(baseline(7e06:7.4e06));
figure; plot(BaselineSample)
[f fftdata]=fftshow(double(BaselineSample));
figure; findpeaks(fftdata, 'MinPeakHeight',1, 'MinPeakDistance',2000); 
% axis([0 20000000 0 10])                                                                                                                                                                                                                                                                                                                                                                                      inf -inf inf])
[pks, locs]=findpeaks(fftdata, 'MinPeakHeight',1, 'MinPeakDistance',2000);
xfreq=f(locs);%xaxis
ypeaks = fftdata(locs);%yaxis

% fftshow'09152016 with pa Fish3' all were below 50Hz, so no filtering done

%% Ripple  Analysis
 
% STEPS
% 1.plot all frequency bands : wideband (JY), filter 500hz, FR, R, Gamma, beta, alpha
% 2.HFO detection:
%         a. Find mean baseline of 10s
%         b. Draw threhold as 3SD(baseline of 10s) )
%         c. detect 4 or more consecutive peaks above threshold
%         d. time lag between peaks 2-4ms=FR, 5-12.5ms=R
%         e. plot R and FR together to check for overlapping event
%         f. exclude overlapping events
%         g. subplot FRs and Rs
% 3. set bin =1s. Identify number of FR per 1s 
% 4. plot them together and average. Do same for R
% 5. plot R and Fr with y axis as rate bin
% 4. spectrogram?
% 6. do for both PTZ and PB
% 7. mesure change in PB amplitude (JY) and frequency wrt PTZ

%STEP 

ephys= double(postPTZ);
% go to fastripplegammaPower(postPTZ)


[c k] = cheby1(2,0.5,[80 200]/5000,'bandpass'); 
A = filtfilt(c,k,ephys); 
threshold= 3*std(A(2.59*10^6:2.69*10^6));
A_phase = angle(hilbert(A));

% -A cause it has more longer peaks on the negative side
figure;findpeaks(-A, 'MinPeakHeight',  threshold, 'MinPeakDistance',5);
hold on
plot(A_phase*100,'r')

[pks locs]=findpeaks(-A, 'MinPeakHeight',  threshold, 'MinPeakDistance',5);
figure; plot(A_phase(locs))
plot(A_phase(locs),'*')
plot(abs(A_phase(locs)),'*')
figure
plot(locs,abs(A_phase(locs)),'*')
hold on
plot(-A/600)% Reducing magnitude by 10x
plot(locs,pks/600,'rs')

% manually find threshold (about >= 2.5) where peaks are correctly distinguishable
hold on 
Thold=2.73
 plot([0 length(A)],[2.73 2.73],'g')

 % excluding locs below threshold
Ax=abs((A_phase(locs)))
locs1=locs

ind=Ax<2.73
Ax(find(ind))=0
locs1(find(ind))=0
hold on
plot(locs1,Ax,'g*')

% excluding zeros and creating new arrays
j= 1;
for i = 1: length(locs)
    if locs1(i)>0 
        locs_N(j) = locs1(i);
        pks_N(j) = pks(i);
        j = j +1;
    end
end
hold on
plot(locs_N,pks_N/600,'ks')

% ripples peaks between 5 to 12.5ms
figure
plot(locs_N,pks_N/600,'gs')
hold on
plot(-A/600) 

B=diff(locs_N')
figure; plot(B)
C=find(B<50 | B>125) 
B(C)=0
plot(B,'r')

j=1
ind=find(B);
for i = 1: length(ind)
       locs_N2(j) = locs_N(ind(i));
        pks_N2(j) = pks_N(ind(i));
        j=j+1;
   
end
hold on
plot(locs_N2,pks_N2/600,'ko')

%finding consecutive peaks in groups of more than 4 
%gap betwen each consecutive spike ~ 100

D= diff(locs_N2')
% figure; plot(D)
 E=find(D>300) 
 D(E)=0
% hold on; plot(D,'r')

j=1
ind=find(D);
for i = 1: length(ind)
       locs_N3(j) = locs_N2(ind(i));
        pks_N3(j) = pks_N2(ind(i));
        j=j+1;
   
end
hold on
plot(locs_N3,pks_N3/600,'r*')



figure; spectrogram(ephys)




figure
subplot(411)
plot(postPTZ(900000:1420000));axis([-inf inf -inf inf])
subplot(412)
plot(P)
subplot(413)
find(f>5,1)
plot(abs(s(5,:)))
title(num2str(f(5)))
subplot(414)
find(f>50,1)
plot(abs(s(46,:)))
title(num2str(f(46)))
h=get(gcf,'children');
linkaxes(h,'x')
subplot(411)
plot(q,postPTZ)
e1 = double(postPTZ(2620*500:2780*500));
fftshow(e1,10000);


figure;  spectrogram(x2,1000, 500, 1024,10000,'yaxis');
[s f t] = spectrogram(x2,1000, 500,4096,10000,'yaxis');
size(s)
size(f)
t(1)
5000/2049
f(1)
f(2)
find(f>60,1)
figure
imagesc(abs(f(1:25,:)))
imagesc(abs(s(1:25,:)))
f(1:25)
[s f t] = spectrogram(x2,1000, 500,9192,10000,'yaxis');
size(s)
[s f t] = spectrogram(x2,500, 250,9192,10000,'yaxis');
size(s)
f(1:50)'
find(f>60,1)
imagesc(abs(s(1:56,:)))
ylabel(f(1:56))
ylable('off')
ylabel('off')
P = sum(abs(s(1:56,:)),1);
figure
plot(p)
plot(P)
 subplot(211)
plot(x2)
subplot(212)
plot(p)
plot(P)
[s f t] = spectrogram(x2,1000, 500,9192,10000,'yaxis');
find(f>60,1)
P = sum(abs(s(1:56,:)),1);
subplot(212)
plot(P_
plot(P)
[s f t] = spectrogram(double(postPTZ),1000, 500,9192,10000,'yaxis');
P = sum(abs(s(1:56,:)),1);
figure
subplot(211)
plot(postPTZ)
subplot(212)
plot(P)
axis([-inf inf -inf inf])
subplot(211)
axis([-inf inf -inf inf])
figure
subplot(411)
plot(postPTZ)
subplot(412)
plot(P)
subplot(413)
find(f>5,1)
plot(abs(s(5,:)))
title(num2str(f(5)))
subplot(414)
find(f>50,1)
plot(abs(s(46,:)))
title(num2str(f(46)))
help linkaxes
h=get(gcf,'children');
linkases(h,'x')
linkaxes(h,'x')
w = downsample(postPTZ,500);
subplot(411)
plot(w)
q = 1:length(postPTZ);
length(q)
size(s)
length(q)/21463
q= q/500;
plot(q,postPTZ)
subplot(414)
plot(sum(abs(s),1))
e = double(postPTZ(2100*500:3250*500));
fftshow(e,10000);
help fftshow
fftshow(e,10000);
e = double(postPTZ(2150*500:2350*500));
fftshow(e,10000);
figure
figure
subplot(411)
plot(postPTZ)
subplot(412)
plot(P)
subplot(413)
find(f>5,1)
plot(abs(s(5,:)))
title(num2str(f(5)))
subplot(414)
find(f>50,1)
plot(abs(s(46,:)))
title(num2str(f(46)))
h=get(gcf,'children');
linkaxes(h,'x')
subplot(411)
plot(q,postPTZ)
e1 = double(postPTZ(2620*500:2780*500));
fftshow(e1,10000);


%%%%%%%%%%%%%%%%%%%%%%%%
figure
[s f t] = spectrogram(double(postPTZ),1000,  500,9192,10000,'yaxis');
P = sum(abs(s(1:56,:)),1);

subplot(411)
q = 1:length(postPTZ);
k=length(q)/length(P); %= 500
q= q/k;
plot(q/20,postPTZ);  xlabel('sec')

subplot(412)
qp = 1:length(P);
plot(qp/20,P);  

subplot(413)
find(f>5,1)
plot(qp/20,abs(s(5,:)));axis([-inf inf -inf inf])
title(num2str(f(5)))
% title('5Hz')

subplot(414)
find(f>50,1)
plot(qp/20,abs(s(46,:))); axis([-inf inf -inf inf])
title(num2str(f(46)))
% plot(qp/20, abs(mean(s(38:46,:)))); axis([-inf inf -inf inf]) % plot 38 to 60~ 40 to 50Hz
%  title('40-50Hz')
% % 
% hold on; 
h=get(gcf,'children');
linkaxes(h,'x')

%%%%%%%%%%%%%%%%%%%%%%%%


postPTZ50Hz=abs(s(46,:));
figure
plot(postPTZ50Hz)

m(1)=1
m(2)=300
% m=ginput(2)
Sz=postPTZ50Hz(m(1)*20:m(2)*20);
Sz1=Sz;
figure; subplot(211)
plot(Sz)
subplot(212)
plot(postPTZ(m(1):m(2)*10000))


% m=ginput(4)
% base=Sz(m(1):m(2));
% seizureBase=Sz(m(3):m(4));
% 
% stdB=std(base)
% stdSB=std(seizureBase)
% % hold on; plot([1,length(Sz)],[ 3*stdB  3*stdB], 'k-')
% hold on; plot([1,length(Sz)],[ 2*stdSB  2*stdSB])

Sz=Sz1;
kk=find(Sz<(2*stdSB));
Sz(kk)=nan;
 
figure
plot(Sz1);
hold on
plot(Sz)


xx=find(isfinite(Sz)');
figure; plot(Sz1(xx(1):xx(length(xx))))

Szduration=(xx(length(xx))-xx(1))/10 % in sec
Szmaxamp=max((Sz1(xx(1):xx(length(xx)))));
 

xx=find(isfinite(Sz)');
xx1=diff(xx)
figure; plot(xx1)
k = [true;xx1>5];
hold on; plot(k)
xx1(k)=nan;
hold on; plot(xx1)

figure;plot(Sz(xx))
hold on; plot(Sz)


A =  Sz;
pos = [true, isnan(A(:, 1)).', true];
ini = strfind(pos, [true, false]);
fin = strfind(pos, [false, true]) - 1;
C   = cell(1, length(ini));
for iC = 1:length(ini)
C{iC} = A(ini(iC):fin(iC), :);
%   C{iC} = A((fin(iC)+1):(ini(iC+1)-1), :);

end

figure;plot(C)


%
 
A = [Nan,NaN;200,NaN;NaN,NaN;-12,14;-13,14;-14,13;-15,13;NaN,NaN;NaN,NaN;...
      NaN,NaN;NaN,NaN;NaN,NaN;-15,13;-14,13;-13,14;-12,14;-12,13;-11,13;...
      -10,14;-9,14;NaN,NaN;NaN,NaN;NaN,NaN;NaN,NaN;-11,12;-13,13;NaN,NaN];
pos = [true, isnan(A(:, 1)).', true];
ini = strfind(pos, [true, false]);
fin = strfind(pos, [false, true]) - 1;
C   = cell(1, length(ini));

%
for i = 1:length(ini)
 if ini(i)==fin(i)
    A(ini(i)) =NaN;
 end
end
pos = [true, isnan(A(:, 1)).', true];
ini = strfind(pos, [true, false]);
fin = strfind(pos, [false, true]) - 1;
C   = cell(1, length(ini));

%
for iC = 1:length(ini)-1
 nn(iC)= ini(iC+1)-fin(iC);
 end
 hold on; plot(nn)

 %
for iC = 2:length(ini)
    if ini(iC+1)-fin(iC)<15
        fin(iC)=fin(iC+1);
        iC=iC-1
    end
end
 


 
 %
for iC = 1:length(ini)
    if ini(iC+1)-fin(iC)<20
     for jC=iC:length(ini)
          if ini(jC+1)-fin(jC)<20
        jC{iC} = A(ini(iC):fin(iC+1), :);
   else
  C{iC} = A(ini(iC):fin(iC), :);
    end
end



 
 %
for iC = 1:length(ini)-1
    if ini(iC+1)-fin(iC)<20
           fin(iC)=fin(iC+1);
        C{iC} = A(ini(iC):fin(iC+1), :);

    else
  C{iC} = A(ini(iC):fin(iC), :);
    end
end




%
sigden = cmddenoise(Sz,'sym5',5,'s',2);
plot(Sz);
hold on;
plot(sigden,'r','linewidth',2);
axis tight;

%
figure; findpeaks(Sz, 'MinPeakHeight',3*stdB)
 
% find 
Sz2(1:5)
for i=1:5:length(Sz)
 if length(find(isnan(Sz2(i:i+4))))>length(find(isfinite(Sz2(i:i+4))))
Sz3(i:i+4)=nan
 else
     Sz3(i:i+4)= Sz2(i:i+4)
 end
end

 
 


 
 
j=1
n=find(Sz);
for i = 1: length(Sz)
    if
    Sz(j)=Sz(i)
    
       locs_N2(j) = locs_N(ind(i));
        pks_N2(j) = pks_N(ind(i));
        j=j+1;
   
end 

%%%%%%%%%%%%%%%%%%%%%%%%
 















x2=double(postPTZ1);
figure;
subplot(414)
imagesc((abs(s(1:56,:))))
imag=(abs(s(1:56,:)));
imgMirror = flipdim(imag,1);
imagesc(imgMirror)
NumTicks = 4;
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks))
set(gca,'YTickLabel',{'60', '40','20', '0'})
ylabel('Hz')

subplot(413)
[s f t] = spectrogram(x2,1000, 500,9192,10000,'yaxis');
find(f>50,1) %=46
plot(abs(mean(s(38:46,:)))); axis([-inf inf -inf inf]) % plot 38 to 60~ 40 to 50Hz
% hold on ; plot(abs(s(46,:))); axis([-inf inf -inf inf])
title('40-50Hz')

title(num2str(f(46)))

subplot(412)
[s f t] = spectrogram(x2,1000, 500,9192,10000,'yaxis');
find(f>5,1)
plot(abs(s(5,:)));axis([-inf inf -inf inf]);title('5Hz')

% title(num2str(f(5)))
% P = sum(abs(s(1:56,:)),1);
% plot(P)

subplot(411)
q=1:length(x2);
plot(q/10000,x2);xlabel('sec')
%%%%%%%%%%%%%%%%%%%%%%%

 
 
 