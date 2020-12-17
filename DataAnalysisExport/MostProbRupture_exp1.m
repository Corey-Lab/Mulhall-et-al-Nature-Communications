function [Frup,dF_all,tau_all,F_all,dtau_all,lr,tau_ASH_all,F_ASH_all,dtau_ASH_all] = MostProbRupture_exp1(F,showp)
%Computes the most probable rupture force from a set of data by binning the
%data data and finding the bin with the highest counts, the bin widths are
%chosen using the rule for fiting guassian distributions. The bin centers
%are then scanned to yield the bin center that has the highest counts.
%   Detailed explanation goes here

if showp==1
    figure;
end


tau_all=[];
tau_ASH_all=[];
dtau_ASH_all=[];
F_ASH_all=[];
F_all=[];
dtau_all=[];

for i=1:length(F)
    
    
    Ft=F{i}.data(:,2);
    %% Method 1 Choose the "best" bin center
    dF=(prctile(Ft,75)-prctile(Ft,25))/length(Ft)^(1/3);
%     dF=3.5*std(Ft)/length(Ft)^(1/3);
    Fmin=min(Ft);
    Fmax=max(Ft);
    Nj=6;
    counts=[];
    centers=[];
    
    for j=1:Nj
        
        [counts,centers]=hist(Ft,[Fmin-dF+j/Nj*dF:dF:Fmax+dF]);
        
        [cj(j),fmax]=max(counts);
        fj(j)=centers(fmax);
        
    end
    
    
    dF_all(i,1)=dF;
    
    [cnt,jmax]=max(cj);
    
    Frup(i,1)=fj(jmax);
    
    [cnts,f]=hist(Ft,[Fmin-dF+jmax/Nj*dF:dF:Fmax+dF]);
    
    %% Method 2 Create Cumulative distribution and fit
    [f_e,x_e]=ecdf(Ft);
    ag=0;
    bg=median(Ft);
    cg=std(Ft);
    fo=fitoptions('Method','NonlinearLeastSquares','MaxIter',800,'MaxFunEvals',1200,...
        'StartPoint', [ag,bg,cg], ...
        'Lower', [0, 0, 0], ...
        'Upper', [1,max(Ft)*2, max(Ft)*10]);
    ft=fittype('a/(1+exp(-(x-b)/c))','options',fo);
    [fd,fs,fo]=fit(x_e,f_e,ft);
    dFit=diff(confint(fd))/2;
    
    Frup(i,2)=fd.b;
    dF_all(i,2)=dFit(2);
    if i==8
        
        10
        
    end
    
    %% Method 3 Use the ksdensity function
    [p_k,f_k,dF_all(i,3)]=ksdensity(Ft);
    [maxp,max_index]=max(p_k);
    Frup(i,3)=f_k(max_index);
    
    
    %% Method 4 What does "best" bin center really mean?? I don't know. Use the average shifted histogram.
    
    [Frup(i,4),dF_all(i,4),ASH] = AverageShiftedHistogram(Ft,4);
    
    %% Create Lifetime data from Force histrogram.
    
    
    lrt=F{i}.lr;
    %     Nmax=sum(cnts);
    %     Nl(1)=Nmax;
    lr(i)=lrt;
    %
    %     for j=2:length(f)
    %         Nl(j)=Nl(j-1)-cnts(j-1);
    %     end
    %     frem=find(cnts==0 | cnts==1);
    %
    %     cnts(frem)=[];
    %     f(frem)=[];
    %     Nl(frem)=[];
    %     p=cnts./Nl;
    %     dp=sqrt(1-p).*p;
    %
    %     tau_f=dF/lrt./p;
    %     dtau_f=dp*dF/lrt./(p.^2);
    %
    %
    %     tau_all=[tau_all;tau_f'];
    %     F_all=[F_all;f'];
    %     dtau_all=[dtau_all;dtau_f'];
    %
    
    
    [F_all,tau_all,dtau_all]=LifetimeConvert(lrt,cnts,f,F_all,tau_all,dtau_all);
    
    %% Create Lifetime data from Average Shifted Histogram
    
    f_ASH=ASH(:,1)';
    cnts_ASH=round(ASH(:,2))';
    
    [F_ASH_all,tau_ASH_all,dtau_ASH_all]=LifetimeConvert(lrt,cnts_ASH,f_ASH,F_ASH_all,tau_ASH_all,dtau_ASH_all);
    
    %%
    
    if showp==1
        
        subplot(5,1,1);
        plot(fj,cj,'bo');
        xlabel('Rupture Force');
        ylabel('Counts');
        subplot(5,1,2);
        hold off;
        stairs(f,cnts/sum(cnts)/(f(2)-f(1)),'r','LineWidth',2);
        hold on;
        stairs(ASH(:,1),ASH(:,2)/sum(ASH(:,2))/(ASH(2,1)-ASH(1,1)),'b','LineWidth',2);
        plot(f_k,p_k);
        syms x
        cdf_fit = fd.a/(1+exp(-(x-fd.b)/fd.c));
        pdf_fit = diff(cdf_fit,x);
        pdf_fun = matlabFunction(pdf_fit);
        plot(f_k,pdf_fun(f_k),'r');
        
        xlabel('Rupture Force');
        ylabel('pdf');
        
        subplot(5,1,3);
        plot(x_e,f_e);
        hold on;
        plot(fd);hold off;
        xlabel('Rupture Force');
        ylabel('CDF')
        subplot(5,1,4);
        plot(f_k,p_k);
        xlabel('Rupture Force');
        ylabel('pdf');
        subplot(5,1,5);
        hold off;
        semilogy(F_all,tau_all,'ro');
       
        hold on;
        
        
        semilogy(F_ASH_all,tau_ASH_all,'bo');
       
        Frup(i,:)
        input('Enter to continue');
        
    end
    
    
end





end

