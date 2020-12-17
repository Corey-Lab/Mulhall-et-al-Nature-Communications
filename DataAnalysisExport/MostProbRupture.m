function [Frup,dF_all,tau_all,F_all,dtau_all] = MostProbRupture(F,lr,showp)
%Computes the most probable rupture force from a set of data by binning the
%data data and finding the bin with the highest counts, the bin widths are
%chosen using the rule for fiting guassian distributions. The bin centers
%are then scanned to yield the bin center that has the highest counts.
%   Detailed explanation goes here

if showp==1
    figure;
end
tau_all=[];
F_all=[];
dtau_all=[];
for i=1:length(F(:,1))
    
    
    Ft=F(i,:);
    
    dF=3.5*std(Ft)/length(Ft)^(1/3);
    Fmin=min(Ft);
    Fmax=max(Ft);
    
    
    for j=1:20
        
        [counts,centers]=hist(Ft,[Fmin-dF+j/20*dF:dF:Fmax+dF]);
        
        [cj(j),fmax]=max(counts);
        fj(j)=centers(fmax);
        
    end
    
    dF_all(i)=dF;
    
    [cnt,jmax]=max(cj);
    Frup(i)=fj(jmax);
    
    [cnts,f]=hist(Ft,[Fmin-dF+jmax/20:dF:Fmax+dF]);
    lrt=lr(i);
    Nmax=sum(cnts);
    Nl(1)=Nmax;
    
    for j=2:length(f)
        Nl(j)=Nl(j-1)-cnts(j-1);
    end
    frem=find(cnts==0 | cnts==1);
    cnts(frem)=[];
    f(frem)=[];
    Nl(frem)=[];
    p=cnts./Nl;
    dp=sqrt(1-p).*p;
    
    tau_f=dF/lrt./p;
    dtau_f=dp*dF/lrt./(p.^2);
    
    
    tau_all=[tau_all;tau_f'];
    F_all=[F_all;f'];
    dtau_all=[dtau_all;dtau_f'];
    
    
    
    
    if showp==1
        
        subplot(1,2,1);
        plot(fj,cj,'bo','MarkerFaceColor','b');
        subplot(1,2,2)
        hist(Ft,[Fmin-dF+jmax/10:dF:Fmax+dF]);
        input('Enter to continue');
    end
    
    
end





end

