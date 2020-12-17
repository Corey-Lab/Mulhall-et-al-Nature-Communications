function [ Frup,trup,lr ] = RuptureDetect( data,disp,Navg )
%Detect Rupture evenet by looking at derivative in force

F=data.k*(data.xt-data.xb);
dF=diff(F);
t=[1:length(F)]*data.dt;

[dFmax,imax]=min(dF);

Frup=mean(F(imax-Navg:imax));
Fhigh_res=max(data.k*(data.xt_q-data.xb_q));
Frup=max([Frup,Fhigh_res]);
trup=imax*data.dt;

t_fit=[1:imax]*data.dt;
p=polyfit(t_fit,F(1:imax),1);

if disp==1
    hold off
    plot(t,F);
    hold on;
    plot(t,dF(end-1),'r')
    plot(trup,Frup,'go','MarkerFaceColor','g');
    plot(t_fit,polyval(p,t_fit),'r','linewidth',2);
end

lr=p(1);

end

