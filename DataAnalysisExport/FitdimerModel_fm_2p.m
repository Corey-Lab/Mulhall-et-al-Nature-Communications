function [fitout] = FitdimerModel_fm_2p(fbN,koffN,C,fc,konN,model,rdimer,Fdimer,W,method,dna)
%Fit dimer model using fminsearch using for fit 2 fit parameters: kon, koff
%and , Input paramters should be C and fc.
%   Detailed explanation goes here

switch(method)
    
    case(1)
        Xsq=@(x)sum((DimerForceNumericalSolution(rdimer,x(3),x(4),x(1)/1000,x(2),konN,model,0.00000000000001,0)-Fdimer).^2.*W);
    case(2)
        Xsq=@(x)(sum((DimerForceNumericalSolution_Median(rdimer,fbN,x(1),C,fc,x(2),model,dna)-Fdimer).^2.*W)+1000*sum(x<0));
end


[fd]=fminsearch(Xsq,[koffN,konN]);



%% calculate error in chisqr
df3=[0.01,1000];
j=1;

%%
% for i=-3:3
%     Xe(j)=Xsq([fd(1)+df*i,fd(2),fd(2]);
%     f(j)=fd(1)+df*i;
%     j=j+1;
% end
% fz=find(f<0);
% f(fz)=[];
% Xe(fz)=[];
% pf=polyfit(f,Xe,2);
% yf=polyval(pf,f);
% Q1=100*mean(abs(yf-Xe)./Xe);
% e1=2*sqrt(2/pf(1));
%%
% 
% j=1;
% 
% for i=-3:3
%     Xe(j)=Xsq([fd(1),fd(2)+df*i]);
%     f(j)=fd(2)+df*i;
%     j=j+1;
%     
% end
% 
% fz=find(f<0);
% f(fz)=[];
% Xe(fz)=[];
% pf=polyfit(f,Xe,2);
% yf=polyval(pf,f);
% Q2=100*mean(abs(yf-Xe)./Xe);
% e2=2*sqrt(2/pf(1));
% 


for k=1:2
    j=1;
    df=df3(k);
    %loop over the 4 different fit parameters
    for i=-5:5
        %Scan the Chisquare for the kth fit paramter
        fdt=fd;
        fdt(k)=fdt(k)+df*i;
        Xe(j)=Xsq(fdt);
        f(j)=fdt(k)+df*i;
        j=j+1;
    end
    
    fz=find(f<0);
    f(fz)=[];
    Xe(fz)=[];
    pf=polyfit(f,Xe,2);
    yf=polyval(pf,f);
    Q(k)=100*mean(abs(yf-Xe)./Xe);
    err(k)=2*sqrt(2/pf(1));
    
end


%% [fbN,koffN,konN]

% fitout.fb=fd(1);
% fitout.dfb=err(1);
fitout.koff=fd(1);
fitout.dkoff=err(1);
fitout.kon=fd(2);
fitout.dkon=err(2);

fitout.Q=Q;
fitout.Xsq=Xsq(fd);
disp('***');
disp(['Xsqr ',num2str(Xsq(fd))]);
disp(['kon(1/M/sec) : ',num2str(fitout.kon),'+/- ',num2str(fitout.dkon),' Q: ', num2str(fitout.Q(2))]);
%disp(['fb(pN) : ',num2str(fitout.fb),'+/- ',num2str(fitout.dfb),' Q: ', num2str(fitout.Q(1))]);
disp(['koff(1/sec) : ',num2str(fitout.koff),'+/- ',num2str(fitout.dkoff),' Q: ', num2str(fitout.Q(1))]);
disp('***');

end

