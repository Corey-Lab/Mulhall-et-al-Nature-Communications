function [fitout] = FitdimerModel_2Barrier(fb1N,koff1N,fb2N,koff2N,rdimer,Fdimer,W)
%Fit dimer model using fminsearch using for fit paramters C, fc,koff,fb
%   Detailed explanation goes here

       



Xsq=@(x)sum((Dimer_2BarrierSolve(rdimer,x(1),x(2),x(3),x(4))-Fdimer).^2.*W+sum(x<0)*1000);
   

[fd]=fminsearch(Xsq,[fb1N,koff1N,fb2N,koff2N]);



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

fitout.fb1=fd(1);
fitout.dfb1=err(1);

fitout.koff1=fd(2);
fitout.dkoff1=err(2);

fitout.fb2=fd(3);
fitout.dfb2=err(3);


fitout.koff2=fd(4);
fitout.dkoff2=err(4);

fitout.Q=Q;
fitout.Xsq=Xsq(fd);
disp('***');
disp(['Xsqr ',num2str(Xsq(fd))]);
disp(['koff1(1/sec) : ',num2str(fitout.koff1),'+/- ',num2str(fitout.dkoff1),' Q: ', num2str(fitout.Q(1))]);
disp(['koff2(1/sec) : ',num2str(fitout.koff2),'+/- ',num2str(fitout.dkoff2),' Q: ', num2str(fitout.Q(4))]);
disp(['fb1(pN) : ',num2str(fitout.fb1),'+/- ',num2str(fitout.dfb1),' Q: ', num2str(fitout.Q(2))]);
disp(['fb2(pN) : ',num2str(fitout.fb1),'+/- ',num2str(fitout.dfb2),' Q: ', num2str(fitout.Q(3))]);
disp('***');

end

