freq=[1000,100,10,1];

amp=50;
subplot(1,2,1);
hold off;
for i=1:numel(freq)
    
    [ Fs,Fg,Fp,t] = MotorAdaptation(0,amp,freq(i),30 );
    hold on;
    subplot(1,2,1)
    plot(t,-Fg);
    
%     Fga(i)=mean(-Fg);
%     subplot(1,2,2);
%     plot(t,-Fs);
    
end

% subplot(1,2,2);
% hold on
% semilogx(freq,Fga,'o');