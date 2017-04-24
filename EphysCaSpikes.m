
% load ephys
[Filename,openmatpath]=uigetfile('*.mat;','Please select raw ephys file to open');
cd(openmatpath)
load(char(Filename))

figure;
plot(ephys);axis([-inf inf -inf inf])
%raw

% filter

% select portion
ephys1=ephys;
% ephys1=ephys(1:8e05);
% [f fftdata]=fftshow(ephys1,10000);suptitle( 'raw')
 
% 1Hz to 55Hz
[c k] = cheby1(2,0.5,[10 50]/5000,'bandpass'); 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tiff IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[filename pathname] = uigetfile('*.tif', 'Select a .tif file to read');
cd(pathname);
InfoImage=imfinfo(filename);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage); % =nume(InfoImage);

%green

green_tiff=imread(filename, 1);
for ii = 4 :3: NumberImages
temp_tiff = imread(filename, ii);
green_tiff=cat(3,green_tiff, temp_tiff);
end

% figure
% imagesc(mean(double(green_tiff),3))

for ii = 1 :  size(green_tiff,3)
a(ii)=mean(mean(mean(green_tiff(:,:,ii))));
end
b=a-mean(a);
% figure; plot(b);title(filename)


% 1.Plotting ephys and calcium together
figure; 
i=(1:length(ephys))/Fs;
ha(1)=subplot(211); plot(i,q);axis([-inf inf -inf inf]); title(filename,'FontSize', 9);xlabel('Time (s)');

% time(s)= frame number *256/500
k=(1:length(b))*.512; %*(length(ephys)/(Fs*length(b)));
ha(2)=subplot(212);plot(k,b); axis([-inf inf -inf inf]);xlabel('Time (s)')

% ha(4)=subplot(413);plot(k,b); axis([-inf inf -inf inf]); title('Ca filtered');xlabel('Time (s)')

linkaxes(ha, 'x');      % Link all axes in x
% tt=suptitle(filename);set(tt, 'FontSize', 8);

