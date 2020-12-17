

%% 
sz1=size(lifetime_avg);
Np=sz1(1);
lm=mean(lifetime_avg);
err=std(lifetime_avg)/sqrt(Np);
sz=size(lifetime_avg);
rb=mode(RB);
drb=std(RB);
drb=reshape(drb,[9,16])';
rb=reshape(rb,[9,16])';
namp=sz(2);
cm=colormap(parula(namp));
%% 
pnum=2;
gamma=0.12;
row=1;
data0=mean(lm(1,1,:));
for i=1:namp
    
    data=reshape(lm(1,i,:),[1,sz(3)]);
    delta=reshape(err(1,i,:),[1,sz(3)]);
    subplot(row,2,pnum);
    errorbar(freq',data,delta,'o-','MarkerFaceColor',cm(i,:),'Color',cm(i,:),'DisplayName',[num2str(round(amp(i)/gamma)),' nm']);hold on;
   xlim([0.1 20000]);
   ylim([0 1.2*data0]);
    %subplot(2,2,2);
    %[ Fs,Fg,Fp,t] = MotorAdaptationNoDrag(0,amp(i),1,30,adaptation);
    %Fmax=max(-Fg);
    %dF=abs(Fth_dim-Fmax);
    %[a,b]=min(dF);
    %fcrit=Rl_th(b)/Fth_dim(b);
   % errorbar(freq'/fcrit,data/max(data),delta/max(data),'o','MarkerFaceColor',cm(i,:),'Color',cm(i,:),'DisplayName',[num2str(amp(i)),' nm']);hold on;
    %subplot(2,2,4);
    %errorbar(freq',data/max(data),delta/max(data),'o','MarkerFaceColor',cm(i,:),'Color',cm(i,:),'DisplayName',[num2str(amp(i)),' nm']);hold on;
end
subplot(row,2,pnum);
set(gca,'XScale','log','FontSize',15);
xlabel('F (Hz)');
ylabel('Lifetime (s) ');

%% Mean Force plot

ampp=30;
amp_i=find(amp==ampp);
lr=ampp*freq;
data=reshape(lm(1,amp_i,:),[1,sz(3)]);
figure
semilogx(lr,data,'o')

%%
figure

for i=1:numel(amp)
    [ Fs,Fg,Fp,t] = MotorAdaptationNoDrag(0,amp(i),1,30,adaptation);
    Famp(i)=mean(-Fg);
    plot(t,-Fg);hold on;
end

tau=DimerAverageLifetime_Numerical(Famp,fbN,koffN,CN,fcN,konN,model);
    
%%

Rl_th=[0.005,0.01,0.1,0.39,0.7,[1:10],[10:10:100],[100:50:3000],[4000:1000:10000],[10000:10000:100000]];
Fth_dim =DimerForceNumericalSolution(Rl_th,fbN,koffN,CN,fcN,konN,model,0.00000000000001,0);

%%

sz=size(Forceout)

for j=1:sz(2)
    
    for i=1:sz(1)
        
        Force=Forceout{i,j}.Force;
        t=Forceout{i,j}.t;
        f=find(t<=0.2);
        f1=f(end);
        ForceMax(i,j)=max(Force(f1:end));
        
        
        [ Fs,Fg,Fp,t] = MotorAdaptationNoDrag(0,amp(j),freq(i),30,adaptation );
        Np=1/freq(i)/(t(2)-t(1));
        ForceAvg(i,j)=mean(-Fg(f1:f1+Np));
   
        
    end
    
end
 F1k=round(ForceMax(15,:)'*10)/10
 F1kavg=round(ForceAvg(15,:)')
  data=reshape(lm,[9,16]);
data=data'