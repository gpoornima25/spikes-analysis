function  [ephysM, ephysF50,ephysF50F2,ephysF100] = LoadFilterEphys(filename)
load(char(filename))
filename

%%%% LPF 100Hz
fs=10000;t=(1:length(ephys))/fs;
lpFilt = designfilt('lowpassiir', 'FilterOrder', 2,...
    'PassbandFrequency', 100, 'PassbandRipple', 1, ...
    'SampleRate', fs, 'DesignMethod', 'cheby1');
% fvtool(lpFilt)
ephysF100 = filtfilt(lpFilt,double(ephys));
% [f1 fftdata1]=fftshow(ephysF100,fs);close;

%%%% 1Hz to 50Hz
[c k] = cheby1(2,0.5,[1 50]/500,'bandpass');
ephysF50 = filtfilt(c,k,double(ephys));

%2 Hz notch filter specs
Fs=10000; t = (0:length(ephysF50)-1)/Fs;
%2Hz notch
d = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',1,'HalfPowerFrequency2',3, ...
    'DesignMethod','butter','SampleRate',Fs);
%     fvtool(d,'Fs',Fs)
ephysF50F2 = filtfilt(d,ephysF50);
% [f2 fftdata2]=fftshow(ephysF50F2,10000);suptitle( ' 2Hz notch filter ')



%%%%
ephysM= ephys-mean(ephys);
ephysF50=ephysF50-mean(ephysF50);
ephysF50F2=ephysF50F2-mean(ephysF50F2);
ephysF100=ephysF100-mean(ephysF100);

figure;
ax1=subplot(411);plot(t,ephysM); axis([-inf inf -inf inf]);xlabel('sec')
title('raw')
ax2=subplot(412);plot(t,ephysF50); axis([-inf inf -.05 .05]);xlabel('sec')
title('1-50Hz BPF')
ax3=subplot(413);plot(t,ephysF50F2); axis([-inf inf -.05 .05]);xlabel('sec')
title('1-50Hz BPF+2Hz notchF')
ax4=subplot(414);plot(t,ephysF100); axis([-inf inf -inf inf]);xlabel('sec')
title('100Hz LPF')
linkaxes([ax2,ax1,ax3,ax4],'x');
tightfig
end