
function dBdt = odefcn(t,B,lr,fb,koff,C,fc,kon,model,kn)

%added in the kn term as an attempt to model the alternative explanation
%given by reviewer this is an additional unbinding pathway that is not
%force dependent.
dBdt=zeros(2,1);

switch (model)
    
    case(0)
        Cg=C;
    case(1)
        Cg=C*exp(-lr*t/fc);
    case(2)
        Cg=C*1/(1+lr*t/fc);
    case(3)
        Cg=C*exp(-(lr*t/fc)^2);
    case(4)
        Cg=C/sqrt(2*(lr*t/fc)^2+1);
    case(5)
        Cg=(kon*C*(1+erf((-lr*t)/fc/sqrt(2)))^2);
        
end



k2=koff*exp(lr*t/fb/2)+kn;
k1=koff*exp(lr*t/fb)+kn;

dBdt(2)= -2*k2*B(2) + Cg*kon*B(1);
dBdt(1)= 2*k2*B(2) -k1*B(1) -Cg*kon*B(1);

end