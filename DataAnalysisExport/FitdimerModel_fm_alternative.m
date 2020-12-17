function [fitout] = FitdimerModel_fm_alternative(fbN,koffN,kn,rdimer,Fdimer,W)
%Fit dimer model using fminsearch using for fit paramters C, fc,koff,fb
%   Detailed explanation goes here

       



Xsq=@(x)sum((DimerForceNumericalSolution(rdimer,x(1),x(2),0,1,0,3,x(3),0.00000000000001,0)-Fdimer).^2.*W+sum(x<0)*1000);
   

[fd]=fminsearch(Xsq,[fbN,koffN,kn]);



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


for k=1:3
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

fitout.fb=fd(1);
fitout.dfb=err(1);
fitout.koff=fd(2);
fitout.dkoff=err(2);
fitout.kn=fd(3);
fitout.dkn=err(3);


fitout.Q=Q;
fitout.Xsq=Xsq(fd);
disp('***');
disp(['Xsqr ',num2str(Xsq(fd))]);
disp(['koff(1/sec) : ',num2str(fitout.koff),'+/- ',num2str(fitout.dkoff),' Q: ', num2str(fitout.Q(2))]);
disp(['kn(1/sec) : ',num2str(fitout.kn),'+/- ',num2str(fitout.dkn),' Q: ', num2str(fitout.Q(3))]);
disp(['fb(pN) : ',num2str(fitout.fb),'+/- ',num2str(fitout.dfb),' Q: ', num2str(fitout.Q(1))]);

disp('***');

end

