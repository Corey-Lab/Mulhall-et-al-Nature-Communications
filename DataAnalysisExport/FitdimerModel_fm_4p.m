function [fitout] = FitdimerModel_fm_4p(fbN,koffN,CN,fcN,konN,model,rdimer,Fdimer,W,method,dna)
%Fit dimer model using fminsearch using for fit paramters C, fc,koff,fb
%   Detailed explanation goes here

switch(method)
    
    case(1)
        
        Xsq=@(x)sum((DimerForceNumericalSolution(rdimer,x(3),x(4),x(1)/1000,x(2),konN,model,0,0.00000000000001,0)-Fdimer).^2.*W+sum(x<0)*1000);
    case(2)
        Xsq=@(x)sum((DimerForceNumericalSolution_Median(rdimer,x(3),x(4),x(1)/1000,x(2),konN,model,dna)-Fdimer).^2.*W+sum(x<0)*1000);
end


[fd]=fminsearch(Xsq,[CN,fcN,fbN,koffN]);



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

%[CN,fcN,fbN,koffN]
df4=fd/100;
j=1;


for k=1:4
    j=1;
    f=[];
    Xe=[];
    df=df4(k);
    %loop over the 4 different fit parameters
    for i=-3:3
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
fitout.C=fd(1);
fitout.dC=err(1);
fitout.fc=fd(2);
fitout.dfc=err(2);
fitout.fb=fd(3);
fitout.dfb=err(3);
fitout.koff=fd(4);
fitout.dkoff=err(4);

fitout.Q=Q;
fitout.Xsq=Xsq(fd);
disp('***');
disp(['Xsqr ',num2str(Xsq(fd))]);
disp(['C(mM) : ',num2str(fitout.C),'+/- ',num2str(fitout.dC),' Q: ', num2str(fitout.Q(1))]);
disp(['fc(pN) : ',num2str(fitout.fc),'+/- ',num2str(fitout.dfc),' Q: ', num2str(fitout.Q(2))]);
disp(['fb(pN) : ',num2str(fitout.fb),'+/- ',num2str(fitout.dfb),' Q: ', num2str(fitout.Q(3))]);
disp(['koff(1/sec) : ',num2str(fitout.koff),'+/- ',num2str(fitout.dkoff),' Q: ', num2str(fitout.Q(4))]);
disp('***');

end

