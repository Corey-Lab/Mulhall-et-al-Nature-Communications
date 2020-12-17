function [fitout] = FitdimerModel_fm_Hybrid_simultaneous(fb_1,koff_1,kon_1,C,fc,fbN,koffN,konN,model,rdimer,Fdimer,W,method,rdimerH,FdimerH,Wh,dna,fixfb)
%Fit dimer model using fminsearch using for fit paramters kon_2,koff_2,fb_2
%assuming we know the kon_1, koff_1,C,fc and fb_1
%   Detailed explanation goes here

switch(method)
    
    case(1)
        Xsq=@(x)sum((DimerForceNumericalSolution(rdimer,x(3),x(4),x(1)/1000,x(2),konN,model,0.00000000000001,0)-Fdimer).^2.*W);
    case(2)
        if fixfb==0
            Xsq=@(x)(sum((DimerForceNumericalSolution_Hybrid_Median(rdimer,fb_1,x(1),koff_1,x(2),C,fc,kon_1,x(3),model,dna)-Fdimer).^2.*W)+sum(x<0)*1000 ...
                +sum((DimerForceNumericalSolution_Median(rdimerH,x(1),x(2),C,fc,x(3),model,dna)-FdimerH).^2.*Wh));
        else
            Xsq=@(x)(sum((DimerForceNumericalSolution_Hybrid_Median(rdimer,fb_1,fb_1,koff_1,x(1),C,fc,kon_1,x(2),model,dna)-Fdimer).^2.*W)+sum(x<0)*1000 ...
                +sum((DimerForceNumericalSolution_Median(rdimerH,fb_1,x(1),C,fc,x(2),model,dna)-FdimerH).^2.*Wh));
        end
end
% DimerForceNumericalSolution_Median(lr_in,fb,koff,C,fc,kon,model,dna)
df3=[0.001,0.01,100];
if fixfb==0
    xg=[fbN,koffN,konN];
else
    xg=[koffN,konN];
    df3=df3(2:end);
end

[fd]=fminsearch(Xsq,xg);



%% calculate error in chisqr


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

df3=[0.001,0.01,100];
j=1;

for k=1:numel(xg)
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


%%


if fixfb==0
    fitout.fb=fd(1);
    fitout.dfb=err(1);
    
else
    fitout.fb=fb_1;
    fitout.dfb=0;
    Q=[0,Q];
end


fitout.koff=fd(2-fixfb);
fitout.dkoff=err(2-fixfb);
fitout.kon=fd(3-fixfb);
fitout.dkon=err(3-fixfb);


fitout.Q=Q;
fitout.Xsq=Xsq(fd);
disp('***');
disp(['Xsqr ',num2str(Xsq(fd))]);
disp(['kon : ',num2str(fitout.kon),'+/- ',num2str(fitout.dkon),' Q: ', num2str(fitout.Q(3))]);

disp(['fb(pN) : ',num2str(fitout.fb),'+/- ',num2str(fitout.dfb),' Q: ', num2str(fitout.Q(1))]);
disp(['koff(1/sec) : ',num2str(fitout.koff),'+/- ',num2str(fitout.dkoff),' Q: ', num2str(fitout.Q(2))]);
disp('***');

end

