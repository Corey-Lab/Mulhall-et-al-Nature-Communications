function [Frup,dF] = MostProbRupture_Exp(FL,dL,Lc,showp)
%Computes the most probable rupture force from a set of data by binning the
%data data and finding the bin with the highest counts, the bin widths are
%chosen using the rule for fiting guassian distributions. The bin centers
%are then scanned to yield the bin center that has the highest counts. Also
%filters out data where Lc-dL<contour length <Lc+dL. Force data should be
%in column 2, and contour length is column 1. showp==1 will plot the
%histogram and bin center scan.




    
    Ft=FL(:,2);
    L=FL(:,1);
    
    i_remove=find(abs(Lc-L)<=dL)
    Ft(i_remove)=[];
    
    dF=3.5*std(Ft)/length(Ft)^(1/3);
    Fmin=min(Ft);
    Fmax=max(Ft);
    
    
    for j=1:10
        
        [counts,centers]=hist(Ft,[Fmin-dF+j/10*dF:dF:Fmax+dF]);
        
        [cj(j),fmax]=max(counts);
        fj(j)=centers(fmax);
    end
   
    
    [cnt,jmax]=max(cj);
    
    Frup(i)=fj(jmax);
    
    if showp==1
        figure;
        subplot(1,2,1);
        plot(fj,cj,'bo','MarkerFaceColor','b');
        subplot(1,2,2)
        hist(Ft,[Fmin-dF+jmax/10:dF:Fmax+dF]);
        input('Enter to continue');
    end
    
    


end

