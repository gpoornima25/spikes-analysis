

%%%%%% read ephys
[filename,openmatpath]=uigetfile('*.mat;','Please select raw ephys file to open');
cd(openmatpath);load(char(filename))
[ephysM, ephysF50,ephysF50F2,ephysF100] = LoadFilterEphys(filename);

%%%%%% read tiff
[filename pathname] = uigetfile('*.tif', 'Select a .tif file to read');
cd(pathname); InfoImage=imfinfo(filename);
% [green_tiff,red_tiff,gTrace,rTrace]  = ReadTIFFfast1(filename);
%[green_tiff,red_tiff,gTrace,rTrace]  = ReadTIFFfast2(filename);
[green_tiff,red_tiff,gTrace,rTrace] = ReadTIFFslow(filename);


%% %%% Plotting ephys and calcium together
h=figure;set(h,'position',[480 535 1190 360]);

subplot(4,2,[1,3]);imagesc(mean(double(green_tiff),3));title('green');colorbar
colormap jet;minC=0;maxC=800;set(gca,'clim',[minC,maxC]);axis off
subplot(4,2,[5,7]); imagesc(mean(double(red_tiff),3));title('red');colorbar
colormap jet;minC=0;maxC=800;set(gca,'clim',[minC,maxC]);axis off

Fs=10000;fi=(1:length(ephysM))/Fs;max(fi)
% bk=(1:length(rTrace))*0.002*32; %  1/(0.002*32)= 16Hz
bk=(1:length(rTrace))*0.002*128;%   =4Hz
ha(1)=subplot(4,2,2);plot(fi,ephysM,'k'); axis([-inf inf -inf inf]);title('LFP');
% set(gca,'Visible','off');ax1=gca; ax1.Title.Visible = 'on'; 
ha(2)=subplot(4,2,4);plot(fi,ephysF50,'k'); axis([-inf inf -inf inf]);title('Filtered LFP');
% set(gca,'Visible','off');ax1=gca; ax1.Title.Visible = 'on'; 
ha(3)=subplot(4,2,6);plot(bk,gTrace,'g'); axis([-inf inf -inf inf]);title('F- Green')
ha(4)=subplot(4,2,8);plot(bk,rTrace,'r'); axis([-inf inf -inf inf]);title('F- red')
% set(gca,'Visible','off');ax1=gca; ax1.Title.Visible = 'on'; ;ax1.XAxis.Visible = 'on';xlabel('sec');
linkaxes(ha,'x');xlabel('sec')
suptitle(filename)

%%
%%test pca

[p,q] = rat(3.9063,.0001);%%Since p/q is frame rate in Hz and tolerance of 0.0001 is acceptable
gre=(resample(gTrace,q,p)); 
size(gre);
gre1=(gre(1:end-1)); size(gre1)
size(downsample(ephysM,Fs))



