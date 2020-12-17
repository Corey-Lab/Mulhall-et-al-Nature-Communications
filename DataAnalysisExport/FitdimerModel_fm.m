function [fitout] = FitdimerModel_fm(fbN,koffN,CN,fcN,konN,model,rdimer,Fdimer,W,method,dna)
%Fit dimer model using fminsearch
%   Detailed explanation goes here

switch(method)
    
    
    case(1)
        %Xsq=@(x)sum((DimerForceNumericalSolution(rdimer,fbN,koffN,x(1)/1000,x(2),konN,model,0,0.00000000000001,0)-Fdimer).^2.*W);
        xg=[CN,fcN];
        Xsq=@(x)sum((DimerForceNumericalSolution(rdimer,fbN,koffN,x(1)/1000,x(2),konN,model,0,0.00000000000001,0)-Fdimer).^2.*W+sum(x<0)*1000);
        
    case(2)
        Xsq=@(x)(sum((DimerForceNumericalSolution_Median(rdimer,fbN,koffN,x(1)/1000,x(2),konN,model,dna)-Fdimer).^2.*W)+sum(x<0)*1000);
        
        xg=[CN,fcN];
        
    case(3)
        %Fix CN
        xg=[CN,fcN];
        Xsq=@(x)sum((DimerForceNumericalSolution(rdimer,fbN,koffN,x(1)/1000,x(2),konN,model,0,0.00000000000001,0)-Fdimer).^2.*W+sum(x<0)*1000 + abs(x(1)-CN)*1000);
end



disp('Fitting data')
options = optimset('Display','iter','TolFun',0.01,'TolX',0.01);
[fd]=fminsearch(Xsq,xg,options);

%% calculate error in chisqr
% df=0.1;
% j=1;
% for i=-10:10
% Xe(j)=Xsq([fd(1)+df*i,fd(2)]);
% f(j)=fd(1)+df*i;
% j=j+1;
% end
% fz=find(f<0);
% f(fz)=[];
% Xe(fz)=[];
% pf=polyfit(f,Xe,2);
% yf=polyval(pf,f);
% Q1=100*mean(abs(yf-Xe)./Xe);
% e1=2*sqrt(2/pf(1));
%
% j=1;
% df=0.1;
% f=[];
% Xe=[];
% for i=-3:3
% Xe(j)=Xsq([fd(1),fd(2)+df*i]);
% f(j)=fd(2)+df*i;
% j=j+1;
% end
%
% fz=find(f<0);
% f(fz)=[];
% Xe(fz)=[];
% pf=polyfit(f,Xe,2);
% yf=polyval(pf,f);
% Q2=100*mean(abs(yf-Xe)./Xe);
% e2=2*sqrt(2/pf(1));

%%
df3=fd/100;
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
fitout.Q=Q;
fitout.Xsq=Xsq([fd(1),fd(2)]);
disp('***');
disp(['Xsqr ',num2str(Xsq([fd(1),fd(2)]))]);
disp(['C(mM) : ',num2str(fitout.C),'+/- ',num2str(fitout.dC),' Q: ', num2str(fitout.Q(1))]);
disp(['fc(pN) : ',num2str(fitout.fc),'+/- ',num2str(fitout.dfc),' Q: ', num2str(fitout.Q(2))]);
disp('***');

end

