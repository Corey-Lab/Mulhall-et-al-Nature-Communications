function [ Fs,Fg,Fp,t] = MotorAdaptationNoDrag(excitation,amp,freq,xps,adaptation )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% 
% B=0.06/50; %pN sec/nm
% ks=0.73; %pN/nm
% kg=0.6;
% S=1/12*1000; %1/sec
% C=1900;%nm/sec
% 
% kp=10;%pN/nm probe spring
% 
% Frest=C/S*kg;
% disp(['Resting tension of ',num2str(Frest),' pN']);
% 
% %xps=10;
% 
% xst=(xps-Frest/kp)/(1+ks/kp);
% xmt=C/S-xst;

%%
if adaptation==1
    
    B=0.06/50; %pN sec/nm
    ks=0.73; %pN/nm
    kg=0.6;
    S=1/12*1000*1.1696; %1/sec
    Cl=1900/1.1696;%nm/sec
    
    kp=10;%pN/nm probe spring
    
    Frest=Cl/S*kg;
    disp(['Resting tension of ',num2str(Frest),' pN']);
    xps=30;
    xst=(xps-Frest/kp)/(1+ks/kp);
    xmt=Cl/S-xst;
    
    %[ Fs,Fg,Fp,t] = MotorAdaptation(0,amp,freq,xps);
    
   % Fmax=max(-Fg);
    %disp(num2str(Fmax));
    
else
    B=0.06/50; %pN sec/nm
    ks=0.73; %pN/nm
    kg=0.6;
    kp=10;
    Cl=0;
    S=0;
    Frest=10;
    
    xst=Frest/kg;
    xps=xst+Frest*(1+ks/kg)/kp;
    xmt=0;
  
    Fmax=amp*ks;
    
end




%% 


Tend=max([8/freq,0.5])*2;


switch(excitation)
    
    case(1)
        dt=B/kp/20;
        t=0:dt:Tend;
        
        %Constant Position
        xp=zeros(numel(t),1)+amp+xps;
    case(0)
        
        %sinusoidal
        dt=min([1/freq/500]);
        t=[0:dt:Tend]';
        xp=amp*sin(2*pi*freq*t)+xps;
        
        
end

xs=zeros(numel(t),1);
xg=zeros(numel(t),1);
for i=1:numel(t)
    
    xpt=xp(i);
    
    [xm_out,xs_out,Fgt,Fst,Fpt] = MotorUpdateNoDrag( Cl,S,B,ks,kp,kg,dt,xmt,xpt,xst );
    
    Fg(i)=Fgt;
    
    xst=xs_out;
    
    Fs(i)=Fst;
    
    Fp(i)=Fpt;
    
    xmt=xm_out;
    
    
    
end



end

