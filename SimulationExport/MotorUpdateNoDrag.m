function [ xm_out,xs_out,Fgt,Fst,Fpt ] = MotorUpdateNoDrag( C,S,B,ks,kp,kg,dt,xmt,xpt,xst )
%Motor Update Updates the positions of the moters and springs in the
%adaptation model
%   Detailed explanation goes here

    
    %Calculate Fs
    
    Fst=-ks*xst;
    %Calculate Fs
    
    Fpt=kp*(xpt-xst);
    
    %Calculate Fg
    
    xgt=xst+xmt;
    
    if xgt>=0
        
        Fgt=-kg*(xgt);
        vm=C-S*xgt;
        dxm=vm*dt;
    else
        Fgt=0;
        vm=C;
        dxm=vm*dt;
        
    end
    
    %dxs=dt*(Fgt+Fst+Fpt)/B;
    
    %update new positions and store data
    
    xst=(kp*xpt-kg*xmt)/(kg+ks+kp);
    
    xs_out=xst;
    
    xm_out=xmt+dxm;
    
    
end

