


function [Sname]   = makeMovie(filename, frames,LFP)

%%%%%% movie with ephys
% frames =maskedgreenImage;
% LFP=ephys;
ln = length(filename);
Sname = [filename(1:ln-4) '__MOVIE']
hb=(1: length(LFP))/10000;
writerObj = VideoWriter(Sname);
writerObj.FrameRate = 512/1000;
writerObj.Quality=100;
open(writerObj);
h = figure;
colormap('default')
set(h,'position',[334 62 859 624]);
for i = 1:length(frames)
    t = (.002*128)*i
    h1 = subplot(311);
    %     set(h1,'position',[30/316 336/506 256/316 100/506]);
    plot(hb,LFP); axis([-inf inf -inf inf]);
    hold on
    plot([t t],[-1 -2],'r');
    hold off;
    title(['Time: ' num2str(t) ' ms']);
    
    h2 = subplot(312);
    %     set(h2,'position',[30/316 30/506 256/316 256/506]);
    imagesc(frames(:,:,i));
    title('green');colormap jet;axis off;
    %   minC=0;maxC=5000;set(gca,'clim',[minC,maxC]);axis off
    
    
    h3=subplot(313)
    %     set(h3,'position',[30/316 30/506 256/316 256/506]);
    imagesc(framesRed(:,:,i));
    title('red');colormap jet;axis off;
    %     minC=0;maxC=5000;set(gca,'clim',[minC,maxC]);axis off
    
    
    frame =getframe(h);
    writeVideo(writerObj,frame);
    pause(0.02);
    clf(h);
end
close(writerObj);

end