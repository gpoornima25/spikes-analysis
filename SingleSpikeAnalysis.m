funtion spike 


cr=.002*32;% 32 pixels/line and 0.02  ?
st=170/cr
ed=200/cr
figure;plot(b(st:ed))
e2=eo;
%%
% I1.  FIND CW, SUBTRACT 
% setting CW baseline to 2sec = 32 datapts
%and SMOOTHING (6 datapts~3sec) 
clear s00 cw newcw
 
%%
for i=1:size(e2(st:ed,:),2)
 if  isnan(max(e2(st:ed,i)))
     s00(:,i)=e2(st:ed,i);
     cw(:,i)=s00(1:32,i);
 else
    cw0(:,i)=e2(st:st+15,i); % for 15 sec window
    s00(:,i)=e2(st:ed,i)-nanmean(cw0(:,i));
    cw(:,i)=s00(1:15,i);
    figure;plot(s00(:,i));title(i)
    h=vline(32,'k');
%     h=vline(1,'k');
%     h=vline(EPt,'r');
%     h=hline(nanmean(newcw(:,i)),'k');
  end
end
%cw0(:,i)=e2(t0(i)-35:t0(i)-25); % 10 sec window 

end
