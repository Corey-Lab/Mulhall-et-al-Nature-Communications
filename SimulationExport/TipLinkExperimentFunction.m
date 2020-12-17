function [lifetime_avg,Forceout,Rebind_avg] = TipLinkExperimentFunction(CN,fcN,model,fbN,koffN,konN,Nreps,adaptation,freq,amp,drag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%freq=[1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,1000,10000];
%freq=100;
%freq=fliplr(freq);
lifetime=[];
Rebind=[];
dF=[];

%CN=0.4799/1000
%fcN=10.4
%model=3;
%fbN=14.3;
%koffN=0.58;
%konN=70000;
FN=0;
%Nreps=5000;
phiN=1;
Frest=[1:10:50];
Nmax=length(freq);
lifetime_avg=[];
lifetime_err=[];
FNl=[1:10:50];
FNl=[0,1,5,[10:10:100]];
Frest=10;
%amp=FNl/0.7;
%adaptation=0;
%amp=fliplr(amp);
%%
Forceout={};
tic
for i=1:numel(freq)
    lifetimeExp=[];
    %     FN=FN_Frest(2,j);
    %     Frest=FN_Frest(1,j);
    tic
    for j=1:numel(amp)
        ampN=amp(j);
        data=[];
        [lifetime_t,escape_t,Fout,tout]=TipLinkSimulation_Sine_Adaptation2(koffN,fbN,konN,CN,fcN,freq(i),phiN,ampN,Nreps,1,model,adaptation,drag);
        lifetimeExp(:,j)=lifetime_t';
        Rebind(:,j)=(escape_t(:,1)-1)+escape_t(:,2)-1;
        data.Force=Fout;
        data.t=tout;
        Forceout{i,j}=data;
        
    end
    toc
    
    lifetime_avg(:,:,i)=(lifetimeExp);
    Rebind_avg(:,:,i)=(Rebind);
    %lifetime_err(i,:)=std(lifetimeExp)/sqrt(Nreps);
    
end

toc

end

