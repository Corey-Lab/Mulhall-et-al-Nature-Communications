function [Frup] = MostProbRupture_ks(F,showp)
%Computes the most probable rupture force from a set of data by binning the
%data data and finding the bin with the highest counts, the bin widths are
%chosen using the rule for fiting guassian distributions. The bin centers
%are then scanned to yield the bin center that has the highest counts.
%   Detailed explanation goes here

if showp==1
    
    figure
    
end
 
for i=1:length(F(:,1))
    
    
    Ft=F(i,:);
    [p,f]=ksdensity(Ft);
    [maxp,max_index]=max(p);
    Frup(i)=f(max_index);
    
    
    
    if showp==1
        
        subplot(2,1,1);
        hist(Ft);
        subplot(2,1,2)
        plot(f,p);
        input('Enter to continue');
    end
    
    
end





end

