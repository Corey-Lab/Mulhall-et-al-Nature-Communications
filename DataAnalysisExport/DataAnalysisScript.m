


%% Enter Dimer Data

%Enter the number of rates you have
N=5;
dimer={};

for i=1:N
    
    data=input('Paste data from excel in []: ');
    lr=input('Enter loading rate [pN/sec]: ');
    
    temp.data=data;
    temp.lr=lr;
    
    dimer{i}=temp;
    
end



%% Calculate Most Probable Rupture Force of Dimer data using different methods
%Column1 [Standard Histogram] Column2 [CDF fit Column3 ][KSDensity] Column4
%Average Shifted Histogram method.

showp=1; % Set this to 1 to show plots of the histograms
[F_dimer,dF_dimer,tau_dimer,Ft_dimer,dtau_dimer,lr_dimer,histdata_dimer] = MostProbRupture_CV(dimer,showp,0);
disp('Dimer Data Processed');
%data columns are different methods of analyzing the data {'Histogram','Cdf','ksdensity','Average Shifted Histogram','Median'};
% the first three are most probable rupture force using different methods,
%Standard histogram, Cumulative distrubution function, ksdensity, Average
%shifted histogram, Median is the median rupture force
%% Enter Monomer data for this condition

%Enter the number of rates you have

N=5;

monomer={};


for i=1:N
    
    data=input('Paste data from excel in []: ');
    lr=input('Enter loading rate [pN/sec]: ');
    
    temp.data=data;
    temp.lr=lr;
    monomer{i}=temp;
    
end


%% Calculate Most Probable Rupture Force of Monomer data using different methods
%Column1 [Standard Histogram] Column2 [CDF fit Column3 ][KSDensity] Column4
%Average Shifted Histogram method.

showp=0;

[F_mon,dF_mon,tau_mon,Ft_mon,dtau_mon,lr_mon,histdata_mon] = MostProbRupture_CV(monomer,showp,0);

disp('Monomer Data Processed');

%% Define Parameters for fitting
param=2;% number of parameters to fit the model
hybrid=0;% 1 if the construct is the hybrid
knN=0;


CN=0.5;% zero forces concentration in mM
fcN=10; %force decay for force dependent concentration in pN
konN=7e4;  % Monomer On Rate 1/s
fbN=14.2; % Monomer Force Scale pN
koffN=0.5; % Monomer Off Rate 1/s
model=3; % model parameter, model 3 is exp(-(F/fcN)^2) model, 0 is constant Concentration

disp('Fit parameters reset');


%% Fit MonomerData First, then dimer data. Fit different Force dependent concnetration models
Frest=[[0:0.5:29],[30:10:100],150];
Rl_th=[0.005,0.009,0.1,0.2,0.3,0.39,0.7,[1:10],[10:10:100],[100:100:1000],[1000:1000:10000]];
Rl_th_mon=[[5:10:100],[100:100:1000],[1000:1000:10000]];
disp([num2str(param),' fit parameters']);
Mondata=1;
dnames={'Histogram','Cdf','ksdensity','Average Shifted Histogram','Median'};
fitdata={};
datatype=3; % Use the ksdensity method 
Name=dnames{datatype};
Fdimer=F_dimer(:,datatype);
wdimer=1./dF_dimer(:,datatype).^2;
edimer=dF_dimer(:,datatype);

%plot data
figure;title(['C0_guess',num2str(CN),'mM, fc_guess ' ,num2str(fcN),'pN']);
subplot(1,2,1);
errorbar(lr_dimer,Fdimer,edimer,'o','DisplayName',[Name,' Dimer'])
hold on;


%lr_fit=[lr_mon;lr_dimer];


if datatype==5
    method=2; % if using median rupture force you the model to fit should be the median rupture force model NOT most probable rupture
else
    method=3; 
end

%first Fit monomer data if doing 2 parameter fit;
if (param==2 || param==5)
    
    Fmon=F_mon(:,datatype);
    wmon=1./dF_mon(:,datatype).^2;
    emon=dF_mon(:,datatype);
    
    if param==2
        [fd_mon,fs_mon,fo_mon]=FitMonomerModel(fbN,koffN,lr_mon,Fmon,wmon,method);
        koffN=fd_mon.koff;
        fbN=fd_mon.fb;
    else
        [fd_mon,fs_mon,fo_mon]=FitMonomerModel_alternative(fbN,koffN,knN,lr_mon,Fmon,wmon,method);
        koffN=fd_mon.koff;
        fbN=fd_mon.fb;
        knN=fd_mon.kn;
    end
    
    if method==1
        if param==2
        Fth_mon = fd_mon.fb*log(Rl_th_mon/fd_mon.fb/fd_mon.koff);
        else
           Fth_mon =  fd_mon.fb*log((-2*fd_mon.fb*fd_mon.kn+Rl_th_mon+sqrt(Rl_th_mon.*(-4*fd_mon.fb*fd_mon.kn+Rl_th_mon)))/2/fd_mon.fb/fd_mon.koff);
        end
        
    else
        Fth_mon = fbN*log(1-log(0.5)*Rl_th_mon/koffN/fbN);
    end
    
    subplot(1,2,1);hold on;
    errorbar(lr_mon,Fmon,emon,'s','DisplayName',' Monomer Data');
    subplot(1,2,1);
    plot(Rl_th_mon,Fth_mon,'DisplayName','Monomer Fit');
    tau_mon=1/fd_mon.koff*exp(-Frest/fd_mon.fb);
    
end




models=[3];

for i=1:numel(models)
    
    
    %FitDimer_Monomer(fbN,koffN,CN,fcN,konN,model,lr_m_d,F_m_d,m_index,W)
    %     [fd,fs,fo] = FitDimerModel(fd_mon.fb,fd_mon.koff,CN,fcN,konN,i-1,lr_dimer,Fdimer,wdimer);
    model_i=models(i);
    if datatype==5
        method=2;
    else
        method=3;
    end
    
    switch(param)
        
        case(2)
            fd=FitdimerModel_fm(fbN,koffN,CN,fcN,konN,model_i,lr_dimer,Fdimer,wdimer,method,0);
            knfinal=0;
        case(4)
            
            fd=FitdimerModel_fm_4p(fbN,koffN,CN,fcN,konN,model_i,lr_dimer,Fdimer,wdimer,method,0);
            fd_mon.fb=fd.fb;
            fd_mon.koff=fd.koff;
            knfinal=0;
        case(3)
            fd=FitdimerModel_fm_3pfb(fbN,koffN,CN,fcN,konN,model_i,lr_dimer,Fdimer,wdimer,method,0);
            fd_mon.fb=fbN;
            fd_mon.koff=fd.koff;
            knfinal=0;
        case(5)
            
            %fit the alternative model
            fd=FitdimerModel_fm_alternative(fbN,koffN,kn,lr_dimer,Fdimer,wdimer);
            
       %fd_mon.koff=fd.koff;
           % knfinal=fd.kn;
            fd.C=0;
            fd.fc=1;
    end
    
    if method==2
        
         Fth_dim =DimerForceNumericalSolution_Median(Rl_th,fd_mon.fb,fd_mon.koff,fd.C/1000,fd.fc,konN,model_i,0);
        %Fth_dim =DimerForceNumericalSolution_Median(Rl_th,fd.fb,fd.koff,fd.C/1000,fd.fc,konN,i-1,0);
    else
          
        Fth_dim =DimerForceNumericalSolution(Rl_th,fd_mon.fb,fd_mon.koff,fd.C/1000,fd.fc,konN,model_i,knfinal,0.00000000000001,0);
        %Fth_dim1=DimerForceNumericalSolution(Rl_th,fd.fb,fd.koff,fd.C/1000,fd.fc,konN,model_i,knfinal,0.00000000000001,0);
        Fth_dim_0 =DimerForceNumericalSolution(Rl_th,fd_mon.fb,fd_mon.koff,0,fd.fc,konN,model_i,knfinal,0.00000000000001,0);
        Fth_dim_const =DimerForceNumericalSolution(Rl_th,fd_mon.fb,fd_mon.koff,fd.C/1000,fd.fc,konN,0,knfinal,0.00000000000001,0);
        
        tau_dim=DimerAverageLifetime_Numerical(Frest,fd_mon.fb,fd_mon.koff,fd.C/1000,fd.fc,konN,model);
        tau_dim_0=DimerAverageLifetime_Numerical(Frest,fd_mon.fb,fd_mon.koff,0,fd.fc,konN,model);
        tau_dim_const=DimerAverageLifetime_Numerical(Frest,fd_mon.fb,fd_mon.koff,fd.C/1000,fd.fc,konN,0);
       
      
    end
    subplot(1,2,1);
    plot(Rl_th,Fth_dim,'DisplayName',['Model: ',num2str(model_i),'Co: ',num2str(fd.C),'mM, fc: ',num2str(fd.fc),'pN Xsq: ',num2str(fd.Xsq)]);
    %plot(Rl_th,Fth_dim1,'DisplayName',['Model: ',num2str(model_i),'Co: ',num2str(fd.C),'mM, fc: ',num2str(fd.fc),'pN Xsq: ',num2str(fd.Xsq)]);
    if param==2
        plot(Rl_th,Fth_dim_0,'DisplayName','Model: Zero Concentration');
        plot(Rl_th,Fth_dim_const,'DisplayName',['Constant Concentration,  ',num2str(model_i),'Co: ',num2str(fd.C)]);
    end
    
    subplot(1,2,2);
    plot(Frest,tau_dim,'DisplayName',['Model: ',num2str(model_i),'Dimer Lifetime'],'LineWidth',2);
    hold on;
    
    if param==2
        plot(Frest,tau_mon,'DisplayName','Monomer Lifetime','LineWidth',2);
        plot(Frest,tau_dim_0,'DisplayName',['Model: ',num2str(model_i),'Dimer Lifetime 0 Concentration'],'LineWidth',2);
        plot(Frest,tau_dim_const,'DisplayName',['Model: ',num2str(model_i),'DimerLifetime Constant Concentration'],'LineWidth',2);
    end
    %     fbN=fd_mon.fb;
    %     fcN=fd.fc;
    %     CN=fd.C/1000;
    
    koffN=fd_mon.koff;
    datafit.fd=fd;
    
    %     datafit.fs=fs;
    %     datafit.fo=fo;
    
    datafit.fd_mon=fd_mon;
    if param==2;
        datafit.fs_mon=fs_mon;
        datafit.fo_mon=fo_mon;
    end
    fitdata{i}=datafit;
    fd;
    display(['Model: ',num2str(model_i)]);
    display([sprintf('\r'),'Zero force lifetime :',num2str((fd.C/1000*konN+3*fd_mon.koff)/2/fd_mon.koff^2),'seconds']);
    
end

subplot(1,2,1);
set(gca,'XScale','log','FontSize',20);
xlabel('Loading Rate(pN/s)');
ylabel('Rupture Force');
%title(['C0_{guess}',num2str(CN),'mM, fc_{guess} ' ,num2str(fcN),'pN']);
subplot(1,2,2);
set(gca,'YScale','log','FontSize',20);
xlabel('Force (pN)');
ylabel('Lifetime (seconds)');


%% Fit the homozygous R113G construct by itself

param=2;
Rl_th=[0.005,0.009,0.1,0.2,0.3,0.39,0.7,[1:10],[10:10:100],[100:20:300]];
Mondata=1;
dnames={'Histogram','Cdf','ksdensity','Average Shifted Histogram','Median'};
fitdata={};
datatype=5;
Name=dnames{datatype};
Fdimer=F_dimer(:,datatype);
wdimer=1./dF_dimer(:,datatype).^2;
edimer=dF_dimer(:,datatype);

model=3;
dna=0;

fb_1=15.2;
C=1.1/1000;
koff_1=0.5;
kon_1=70000;
fc=7.97;

method=2;

koffN=1.5;
fbN=15.2;
konN=70000;

switch(param)
    case(3)
        
        
        fd=FitdimerModel_fm_3p(fbN,koffN,C,fc,konN,model,lr_dimer,Fdimer,wdimer,method,dna);
        
    case(2)
        fd=FitdimerModel_fm_2p(fbN,koffN,C,fc,konN,model,lr_dimer,Fdimer,wdimer,method,dna);
        fd.fb=fbN;
end

Ffit=DimerForceNumericalSolution_Median(Rl_th,fd.fb,fd.koff,C,fc,fd.kon,model,dna);

Fpred=DimerForceNumericalSolution_Hybrid_Median(Rl_th,fb_1,fd.fb,koff_1,fd.koff,C,fc,kon_1,fd.kon,model,dna);




%% Fit the hybrid R113G construct

figure;
Rl_th=[0.005,0.009,0.1,0.2,0.3,0.39,0.7,[1:10],[10:5:100],[100:20:300]];
datatype=3;
Fdimer=F_dimer(:,datatype);
wdimer=1./dF_dimer(:,datatype).^2;
edimer=dF_dimer(:,datatype);

simult=1; % simult=1 simultaneiously fits the homozygous and heterozygous data

param=3;
Mondata=1;
dnames={'Histogram','Cdf','ksdensity','Average Shifted Histogram','Median'};
fitdata={};

Name=dnames{datatype};
FdimerH=F_dimerH(:,datatype);% this should be the Homozygous rupture force data.
wdimerH=1./dF_dimerH(:,datatype).^2;
edimerH=dF_dimerH(:,datatype);

subplot(1,2,1);

errorbar(lr_dimer,Fdimer,edimer,'o');hold on;
errorbar(lr_dimerH,FdimerH,edimerH,'s');

model=3;
dna=0;
fb_1=13.5;
C=0.5/1000;
koff_1=0.52;
kon_1=70000;
fc=10.4;
method=2;
koffN=1;
fbN=1.8;
konN=10;
fixfb=0;
fixfc=0;
fcN=10.5;
CN=0.5;
usemedian=0;
fixC=0;

switch(param)
    
    case(3)
        
        switch(simult)
            case(0)
                fd=FitdimerModel_fm_Hybrid(fb_1,koff_1,kon_1,C,fc,fbN,koffN,konN,model,lr_dimer,Fdimer,wdimer,method,dna);
            case(1)
                fd=FitdimerModel_fm_Hybrid__FC_simultaneous(fb_1,koff_1,kon_1,C,fc,fcN,fbN,koffN,konN,CN,model,lr_dimer,Fdimer,wdimer,usemedian,lr_dimerH,FdimerH,wdimerH,dna,fixfb,fixfc,fixC);
        end
        
    case(2)
        fd=FitdimerModel_fm_Hybrid_2p(fb_1,koff_1,kon_1,C,fc,fbN,koffN,konN,model,lr_dimer,Fdimer,wdimer,method,dna);
        fd.fb=fb_1;
end


switch(fixfc)
    
    %DimerForceNumericalSolution_Hybrid_fc(lr_in,fb_1,fb_2,koff_1,koff_2,C,fc_1,fc_2,kon_1,kon_2,model,dna,use_median)
    case(1)
        Ffit=DimerForceNumericalSolution_Hybrid_fc(Rl_th,fb_1,fd.fb,koff_1,fd.koff,C,fc,fd.fc,kon_1,fd.kon*10000,model,dna,usemedian);
        Fpred=DimerForceNumericalSolution(Rl_th,fd.fb,fd.koff,C,fc,fd.kon*10000,model,0,0);
    case(0)
        
        switch(fixC)
            
            case(0)
                Ffit=DimerForceNumericalSolution_Hybrid_fc(Rl_th,fb_1,fd.fb,koff_1,fd.koff,C,fd.C/1000,fc,fd.fc,kon_1,fd.kon*10000,model,dna,usemedian);
                Fpred=DimerForceNumericalSolution(Rl_th,fd.fb,fd.koff,fd.C/1000,fd.fc,fd.kon*10000,model,0,0);
                
            case(1)
                Ffit=DimerForceNumericalSolution_Hybrid_fc(Rl_th,fb_1,fd.fb,koff_1,fd.koff,C,C,fc,fd.fc,kon_1,fd.kon*10000,model,dna,usemedian);
                Fpred=DimerForceNumericalSolution(Rl_th,fd.fb,fd.koff,C,fd.fc,fd.kon*10000,model,0,0);
        end
        
end

plot(Rl_th,Ffit);
plot(Rl_th,Fpred);


subplot(1,2,1);
set(gca,'XScale','log','FontSize',20);
xlabel('Loading Rate(pN/s)');
ylabel('Rupture Force');
%title(['C0_{guess}',num2str(CN),'mM, fc_{guess} ' ,num2str(fcN),'pN']);
subplot(1,2,2);
set(gca,'YScale','log','FontSize',20);
xlabel('Force (pN)');
ylabel('Lifetime (seconds)');

