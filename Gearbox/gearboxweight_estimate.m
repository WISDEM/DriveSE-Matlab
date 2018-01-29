function [Gearboxweight]=gearboxweight_estimate(GearRatio,RotorTorque,Np,GearConfig,Type_ratio,Type_sh)
% subroutine to find speed ratios
% System engineering gearbox sizing tool
% Yi Guo
% NREL NTWC
% 11/9/2012

        if strcmp(GearConfig,'p') ; 
            stagetype=[1];
        elseif strcmp(GearConfig,'e'); 
            stagetype=[2];
        elseif strcmp(GearConfig,'pp');
            stagetype=[1;1];
        elseif strcmp(GearConfig,'ep'); 
            stagetype=[2;1];
        elseif strcmp(GearConfig,'ee');
            stagetype=[2;2];
        elseif strcmp(GearConfig,'eep');
            stagetype=[2;2;1];
        elseif strcmp(GearConfig,'eep_3');
            stagetype=[2;2;1];
        elseif strcmp(GearConfig,'eep_2');
            stagetype=[2;2;1];
        elseif strcmp(GearConfig,'epp');
            stagetype=[2;1;1];
        elseif strcmp(GearConfig,'eee');
            stagetype=[2;2;2];
        elseif strcmp(GearConfig,'ppp');
            stagetype=[1;1;1];
        else
        end
 
[stageratio]=stage_ratio(GearRatio,Np,Type_ratio,GearConfig); % get the speed ratio per stage

N_stages=length(stageratio);        
stagemass=ones(N_stages,1);
stagetorque=ones(N_stages,1);

Ka=0.6; % application factor for weight estimate;
if      (RotorTorque < 2e5);      
            KFact = 850.0;  %KFactor for pitting analysis
elseif (RotorTorque < 7e5)
            KFact = 950.0;
else
            KFact = 1100.0;
end
Kunit=8.029; % unit conversion from Nm to inlb and then lb to N;
  if strcmp(Type_sh,'normal') ; 
      k_sh=1;
  elseif strcmp(Type_sh,'short') ; 
      k_sh=1.25; % taking into account on the shaft length effects on gearbox weight
  end
 %
 stagetorque_0=RotorTorque;
 for i=1:N_stages
 stagetorque(i)=stagetorque_0/stageratio(i);
 stagemass(i)=Kunit*Ka/KFact*stagetorque(i)*stage_mass(stageratio(i),Np(i),stagetype(i));
 stagetorque_0=stagetorque(i);
 end
 Gearboxweight=sum(stagemass)*k_sh;
 
 
 

 
        


 