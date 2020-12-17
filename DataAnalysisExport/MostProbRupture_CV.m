function [Frup,dF_all,tau_all,F_all,dtau_all,lr,histdata] = MostProbRupture_CV(F,showp,hmethod)
%Computes the most probable rupture force from a set of data by binning the
%data data and finding the bin with the highest counts, the bin widths are
%chosen using the rule for fiting guassian distributions. The bin centers
%are then scanned to yield the bin center that has the highest counts.
%   Detailed explanation goes here

if showp==1
    fig1=figure;
end


tau_all=[];
tau_ASH_all=[];
dtau_ASH_all=[];
F_ASH_all=[];
F_all=[];
dtau_all=[];

for i=1:length(F)
    
    sz=size(F{i}.data);
    
    if sz(2)==2
        Ft=F{i}.data(:,2);
    else
        Ft=F{i}.data(:,1);
    end
    %% Method 1 Choose the best bin center/width using Cross Validation
    % From: Nonparametric density estimation and optimal bandwidth selection for protein unfolding and unbinding data.
    dF=1*(prctile(Ft,75)-prctile(Ft,25))/length(Ft)^(1/3);
    %   dF=3.5*std(Ft)/length(Ft)^(1/3);
    Fmin=min(Ft);
    Fmax=max(Ft);
%     if hmethod==1
%         [h_opt,cnts,f] = HistCV(Ft);
%     else
        [h_opt,cnts,f] = HistStd(Ft);
%     end
    [~,jmax]=max(cnts);
    
%     ffit=f;
%     cfit=cnts;
%     
%     fr=find(cfit<cfit(jmax)*0.2);
%     ffit(fr)=[];
%     cfit(fr)=[];
%     
%     if numel(ffit)<3
%         ffit=f;
%         cfit=cnts;
%     end
%     if i==3
%         
%         10
%     end
%     fd_hist=fit(ffit',cfit','gauss1');
%     
%     Frup(i,1)=fd_hist.b1;
     dF_all(i,1)=h_opt/2;
     Frup(i,1)=f(jmax);
   

    histdata.hist=F{i}.data;
    %% Method 2 Create Cumulative distribution and fit
    [f_e,x_e]=ecdf(Ft);
    ag=1;
    bg=median(Ft);
    cg=std(Ft);
    fo=fitoptions('Method','NonlinearLeastSquares','MaxIter',800,'MaxFunEvals',1200,...
        'StartPoint', [ag,bg,cg], ...
        'Lower', [1, 0, 0], ...
        'Upper', [1,max(Ft)*2, max(Ft)*10]);
    ft=fittype('a/(1+exp(-(x-b)/c))','options',fo);
    [fd,fs,fo]=fit(x_e,f_e,ft);
    dFit=diff(confint(fd))/2;
    
    Frup(i,2)=fd.b;
    
    dF_all(i,2)=dFit(2);
    
    
    
    %% Method 3 Use the ksdensity function
    b=KernelCV(Ft,4);
    if hmethod==1
        [p_k,f_k,dF_all(i,3)]=ksdensity(Ft,'Bandwidth',b);
    else
        [p_k,f_k,b]=ksdensity(Ft);
    end
    [maxp,max_index]=max(p_k);
    Frup(i,3)=f_k(max_index);
    dF_all(i,3)=b/2;
    histdata.ks=[p_k',f_k'];
    
    %% Method 4 What does "best" bin center really mean?? I don't know. Use the average shifted histogram.
    
    [Frup(i,4),dF_all(i,4),ASH] = AverageShiftedHistogram(Ft,4);
    histdata.ash=ASH;
    
    
    
    %% Mehtod 5 Use the meidan and median absolute deviation
    Frup(i,5)=median(Ft);
    dF_all(i,5)=1.25*std(Ft)/sqrt(numel(Ft));
    
    %% Create Lifetime data from Force histrogram.
    
    
    lrt=F{i}.lr;
    %     Nmax=sum(cnts);
    %     Nl(1)=Nmax;
    lr(i,1)=lrt;
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
        figure(fig1)
        %subplot(5,1,1);
        % plot(fj,cj,'bo');
        %         xlabel('Rupture Force');
        %         ylabel('Counts');
        subplot(4,1,1);
        hold off;
        bar(f,cnts/sum(cnts)/(f(2)-f(1)),'r','DisplayName','Histogram');
        %stairs(f,cnts,'r','LineWidth',2);
        hold on;
        bar(ASH(:,1),ASH(:,2)/sum(ASH(:,2))/(ASH(2,1)-ASH(1,1)),'b','DisplayName','AverageShiftedHistogram');
        %stairs(ASH(:,1),ASH(:,2),'b','LineWidth',2);
        plot(f_k,p_k,'DisplayName','ksdensity','LineWidth',2);
        %plot(f_k,feval(fd_hist,f_k),'r','LineWidth','2.5');
        syms x
        cdf_fit = fd.a/(1+exp(-(x-fd.b)/fd.c));
        pdf_fit = diff(cdf_fit,x);
        pdf_fun = matlabFunction(pdf_fit);
        %plot(f_k,pdf_fun(f_k),'r');
        
        xlabel('Rupture Force');
        ylabel('pdf');
        
        subplot(4,1,2);
        plot(x_e,f_e);
        hold on;
        plot(fd);hold off;
        xlabel('Rupture Force');
        ylabel('CDF')
        subplot(4,1,3);
        hold off;
        %plot(f,cnts,'bo');
       
        %plot(fd_hist);
         plot(f_k,p_k);
        xlabel('Rupture Force');
        ylabel('pdf');
        subplot(4,1,4);
        hold off;
        semilogy(F_all,tau_all,'ro');
        xlabel('Force [pN]');
        ylabel('liftime [seconds]');
        
        %         hold on;
        %
        %
        %semilogy(F_ASH_all,tau_ASH_all,'bo');
        
        Frup(i,:)
        input('Enter to continue');
        
    end
    
    
end

if showp==1
    set(gca,'FontSize',20);
    
end



end

