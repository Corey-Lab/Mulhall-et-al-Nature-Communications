function [fitout] = FitdimerModel_fm_3pfb(fbN,koffN,CN,fcN,konN,model,rdimer,Fdimer,W,method,dna)
%Fit dimer model using fminsearch using for fit paramters C, fc,koff. fb is
%fixed
%   Detailed explanation goes here

switch(method)
    
    case(1)
        Xsq=@(x)sum((DimerForceNumericalSolution(rdimer,fbN,x(3),x(1)/1000,x(2),konN,model,0.00000000000001,0)-Fdimer).^2.*W);
    case(2)
        Xsq=@(x)sum((DimerForceNumericalSolution_Median(rdimer,fbN,x(3),x(1)/1000,x(2),konN,model,dna)-Fdimer).^2.*W);
end
xg=[CN,fcN,koffN];

[fd]=fminsearch(Xsq,xg);



%% calculate error in chisqr
df=0.01;
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
% 
% 
% for k=1:3
%     j=1;
%     f=[];
%     Xe=[];
%     %loop over the 4 different fit parameters
%     for i=-3:3
%         %Scan the Chisquare for the kth fit paramter
%         fdt=fd;
%         fdt(k)=fdt(k)+df*i;
%         Xe(j)=Xsq(fdt);
%         f(j)=fdt(k)+df*i;
%         j=j+1;
%     end
%     
%     fz=find(f<0);
%     f(fz)=[];
%     Xe(fz)=[];
%     pf=polyfit(f,Xe,2);
%     yf=polyval(pf,f);
%     Q(k)=100*mean(abs(yf-Xe)./Xe);
%     err(k)=2*sqrt(2/pf(1));
%     
% end
% 

%%
df3=fd/100; % start with small guess than increase if needed
j=1;
for k=1:numel(xg)
    j=1;
    df=df3(k);
    pts=5;
    f=zeros(pts*2+1,1);
    Xe=zeros(pts*2+1,1);
    check=0
    
    while check==0
        %loop over the 4 different fit parameters
        for i=1:pts*2+1
            %Scan the Chisquare for the kth fit paramter
            di=i-pts-1;
            fdt=fd;
            ft=fdt(k)+df*di;
            fdt(k)=ft;
            
            Xt= Xsq(fdt);
            Xe(i)=Xt;
            f(i)=ft;
            
        end
        
        fz=find(f<0);
        f(fz)=[];
        Xe(fz)=[];
        Xe=smooth(Xe,3);
        pf=polyfit(f,Xe,2);
        yf=polyval(pf,f);
        
        
        if pf(1)>0
            err(k)=2*sqrt(2/pf(1));
            Q(k)=100*mean(abs(yf-Xe)./Xe);
            check=1;
        else
            df=df*10;
        end
    end
end


%%
fitout.C=fd(1);
fitout.dC=err(1);
fitout.fc=fd(2);
fitout.dfc=err(2);
% fitout.fb=fd(3);
% fitout.dfb=err(3);
fitout.koff=fd(3);
fitout.dkoff=err(3);

fitout.Q=Q;
fitout.Xsq=Xsq(fd);
disp('***');
disp(['Xsqr ',num2str(Xsq(fd))]);
disp(['C(mM) : ',num2str(fitout.C),'+/- ',num2str(fitout.dC),' Q: ', num2str(fitout.Q(1))]);
disp(['fc(pN) : ',num2str(fitout.fc),'+/- ',num2str(fitout.dfc),' Q: ', num2str(fitout.Q(2))]);
%disp(['fb(pN) : ',num2str(fitout.fb),'+/- ',num2str(fitout.dfb),' Q: ', num2str(fitout.Q(3))]);
disp(['koff(1/sec) : ',num2str(fitout.koff),'+/- ',num2str(fitout.dkoff),' Q: ', num2str(fitout.Q(3))]);
disp('***');

end

