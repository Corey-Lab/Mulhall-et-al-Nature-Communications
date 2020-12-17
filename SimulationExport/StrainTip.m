function [] = StrainTip(freq,mag,Frest,k,vmotor)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%lp=2000;


xstart=Frest/k;
xend=0;

dt=1/freq/200;
Tmax=1/freq*2;

t=0:dt:3*Tmax;

xwave=sin(2*pi*freq*t)*mag+xstart;

xendt=zeros(numel(t),1);
for i=1:numel(t)
   
    Ft=k*(xwave(i)-xend);
    
    Ft=Ft*(Ft>0);
    Fwave(i)=Ft;
    if Ft >10
        
        xend=xend+dt*vmotor;
        
    end
    
    if Ft<10
        
        xend=xend-dt*vmotor;
        
    end
    
    xendt(i)=xend;
    
end


figure(1);


subplot(2,1,1);
plot(t,xwave);
hold on;plot(t,xendt,'r');
subplot(2,1,2);
hold on;
plot(t,Fwave);
end

