function [Tbreak,escape_all,Fout,tout] = TipLinkSimulation_Sine_Adaptation2(koff,fb,kon,C,fc,freq,phi,amp,N,Dimer,Cmodel,adaptation,drag)
%Simple Kinetic Simulation involving force dependent off rate, I think this
%is what Welsely was simulating
%Dimer==1 is dimer simulation
%Dimer==0 is monomer simulation
%   Detailed explanation goes here


%% Setup Adaptation parameters

if adaptation==1
    
    B=0.06/50; %pN sec/nm
    ks=0.73; %pN/nm
    kg=0.6;
    S=1/12*1000*1.1696; %1/sec
    Cl=1900/1.1696;%nm/sec
%     
%      S=1/12*1000; %1/sec
%     Cl=2200;%nm/sec
%     
    kp=10;%pN/nm probe spring
    
    Frest=Cl/S*kg;
    disp(['Resting tension of ',num2str(Frest),' pN']);
    xps=30;
    
    
    xst=(xps-Frest/kp)/(1+ks/kp);
    
    xmt=Cl/S-xst;
    
    
    
    if drag==1
        [ Fs,Fg,Fp,tt] = MotorAdaptation(0,amp,freq,xps,adaptation);
    else
        [ Fs,Fg,Fp,tt] = MotorAdaptationNoDrag(0,amp,freq,xps,adaptation);
    end
    
    Fmax=max(-Fg);
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
toff=1/koff;
ton=1/(kon*C);
toff_min=toff*exp(-(Frest+Fmax)/fb);
%shortest time scale we should measure is dependent on maximum force we
%will exhert.


%set simulation shortest time scale in simulation / 50

if freq>0
    if drag==1
    dt=min([toff_min,ton,1/freq,B/kp])/20;
    else
        dt=min([toff_min,ton,1/freq])/20;
    end
    
else
    
    toff_min=toff*exp(-(Frest+Fmax)/fb/2);
    dt=min([toff_min,ton])/10;
    
end



%%



escape_all=zeros(N,2);

parfor i=1:N
    
    state_c1=1;
    state_c2=1;
    bound=state_c1+state_c2;
    escape_c1=0;
    escape_c2=0;
    
    %initialize positions
    if adaptation==1
        xps=30;
        Frest=Cl/S*kg;
        xst=(xps-Frest/kp)/(1+ks/kp);
        xmt=Cl/S-xst;
        
    else
        Frest=10;
        
        xst=Frest/kg;
        xps=xst+Frest*(1+ks/kg)/kp;
        xmt=0;
        
        
    end
    
    
    if Dimer==0
        
        state_c2=0;
        bound=1;
        
    end
    j=1;
    ton=0;
    t=0;
    
    if phi==1
        phi_N=rand*2*pi;
        phi_t(i)=phi_N;
    else
        phi_N=0;
        phi_t(i)=0;
    end
    
    toff=1/koff;
    ton=1/(kon*C);
    toff_min=toff*exp(-(Frest+Fmax)/fb);
    
    
    if freq>0
        
        dt=min([toff_min,ton,1/freq,B/kp])/20;
        
    else
        
        toff_min=toff*exp(-(Frest+Fmax)/fb/2);
        dt=min([toff_min,ton])/10;
        
    end
    
    
    
    
    while bound>=1
        
        t=j*dt;
        %% Update Force and Force dependent rates
        
        %Update Probe position
        xpt=amp*sin(2*pi*freq*t+phi)+xps;
        
        %Update motor/spring positions
        if drag==1
        [ xmt,xst,Fgt,Fst,Fpt ] = MotorUpdate( Cl,S,B,ks,kp,kg,dt,xmt,xpt,xst );
        else
             [ xmt,xst,Fgt,Fst,Fpt ] = MotorUpdateNoDrag( Cl,S,B,ks,kp,kg,dt,xmt,xpt,xst );
        end
        
        F=-Fgt;
        tf=toff*exp(-F/fb/bound);
        
        r1=rand;
        
        prob_off=dt/tf;
        switch (Cmodel)
            
            case(0)
                ton=1/kon/C;
            case(1)
                ton=1/(kon*C*exp(-F/fc));
            case(2)
                ton=1/(kon*C*(1/(1+F/fc)));
            case(3)
                ton=1/(kon*C*exp(-(F/fc)^2));
                
            case(4)
                ton=1/(kon*C/sqrt(2*(F/fc)^2+1));
        end
        prob_on=dt/ton;
        
        %% Evolve State of C1
        
        if state_c1==1
            state_c1=r1>prob_off; % particle becomes unbound if r1<prob
            if state_c1==0
                
                escape_c1=escape_c1+1;
            end
        else
            
            state_c1=r1<prob_on;
        end
        
        %%
        
        
        %% Evolve state of C2 on if Dimer==1 Otherwise this is a monomer simulation
        if Dimer==1
            
            r2=rand;
            if state_c2==1
                state_c2=r2>prob_off;
                if state_c2==0
                    
                    escape_c2=escape_c2+1;
                end
            else
                state_c2=r2<prob_on;
            end
        end
        
        %%
        
        bound=state_c1+state_c2;
        j=j+1;
        if freq==0
            toff_min=toff*exp(-(Frest)/fb/bound);
            dt=min([toff_min,ton])/10;
        end
        
    end
    
    Tbreak(i)=t;
    escape_all(i,:)=[escape_c1,escape_c2];
    
end


%%

 %initialize positions
    if adaptation==1
        xps=30;
        Frest=Cl/S*kg;
        xst=(xps-Frest/kp)/(1+ks/kp);
        xmt=Cl/S-xst;
        
    else
        Frest=10;
        
        xst=Frest/kg;
        xps=xst+Frest*(1+ks/kg)/kp;
        xmt=0;
        
        
    end

Nmax=round(mean(Tbreak)/dt);
Fout=zeros(Nmax,1);
for i=1:Nmax
    t=i*dt;
     %Update Probe position
        xpt=amp*sin(2*pi*freq*t+phi)+xps;
    if drag==1
    [ xmt,xst,Fgt,Fst,Fpt ] = MotorUpdate( Cl,S,B,ks,kp,kg,dt,xmt,xpt,xst );
    else
        [ xmt,xst,Fgt,Fst,Fpt ] = MotorUpdateNoDrag( Cl,S,B,ks,kp,kg,dt,xmt,xpt,xst );
    end
    Fout(i)=-Fgt;
end
tout=dt*[1:Nmax]';

mean(phi_t);
std(phi_t);
end
