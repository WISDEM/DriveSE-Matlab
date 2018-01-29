 function [stageratio]=stage_ratio(GearRatio,Np,type,GearConfig)

% System engineering gearbox sizing tool
% Yi Guo
% NREL NTWC
% 11/9/2012
 %type=1 parallel,type=2, planetary gear;

 
 Kr=0;
 tol=1e-2; % iteration tolerance
 format short e;
if strcmp(type,'empirical')   % Sunderland model to estimate speed ratio
      if strcmp(GearConfig,'p') ; 
            stageratio=[GearRatio];
      elseif strcmp(GearConfig,'e'); 
            stageratio=[GearRatio];
      elseif strcmp(GearConfig,'pp');
            stageratio=[GearRatio^0.5;GearRatio^0.5];
      elseif strcmp(GearConfig,'ep'); 
            stageratio=[GearRatio/2.5;2.5];
      elseif strcmp(GearConfig,'ee');
            stageratio=[GearRatio^0.5;GearRatio^0.5];
      elseif strcmp(GearConfig,'eep');
            stageratio=[(GearRatio/3)^0.5;(GearRatio/3)^0.5;3];
      elseif strcmp(GearConfig,'epp');
            stageratio=[(GearRatio)^(1/3);(GearRatio)^(1/3);(GearRatio)^(1/3)];
      elseif strcmp(GearConfig,'eee');
            stageratio=[(GearRatio)^(1/3);(GearRatio)^(1/3);(GearRatio)^(1/3)];
      elseif strcmp(GearConfig,'ppp');
            stageratio=[(GearRatio)^(1/3);(GearRatio)^(1/3);(GearRatio)^(1/3)];
      else
      end
elseif strcmp(type,'optimal') % iteration based speed ratio calculator
    if strcmp(GearConfig,'epp');
        M1=12;U2_new=1;U2=0;U10=5;U20=4;  %Initial guesses for speed ratios. M1=U1*U2;
             while abs(U2_new-U2)>tol
        par=[M1,Kr,Np(1)];
        U1= fsolve_TE3(@fun_u1_epp,U10,par,tol,1000);  
        U2=M1/U1;
        U3=GearRatio/M1;
        M2=U2*U3;
        par=[M2,Kr,Np(1),GearRatio];
        U2_new= fsolve_TE3(@fun_u2_epp,U20,par,tol,1000);
        M1=M1+0.05;
            end
            stageratio=[U1;U2;U3] ;  
    elseif strcmp(GearConfig,'eep');
         M1=20;U2_new=1;U2=0;U10=5;U20=7;    %Initial guesses for speed ratios. M1=U1*U2;
             while abs(U2_new-U2)>tol
        Kr_1=0;
        Kr_2=0; %2nd stage structure weight coefficient
        par=[M1,Kr_1,Kr_2,Np(1:2),GearRatio];
        U1= fsolve_TE3(@fun_u1_eep,U10,par,tol,1000);  
        U2=M1/U1;
        U3=GearRatio/M1;
        M2=U2*U3;
        par=[M2,Kr_1,Kr_2,Np(1:2),GearRatio];
        U2_new= fsolve_TE3(@fun_u2_eep,U20,par,tol,1000);
        M1=M1+0.05;
            end
            stageratio=[U1;U2;U3];
    elseif strcmp(GearConfig,'eep_3');  % speed ratio calculator when last stage ratio equals three
        Kr_1=0;
        Kr_2=0.8; %2nd stage structure weight coefficient
        M1=GearRatio/3;U10=5;
        par=[M1,Kr_1,Kr_2,Np(1:2),GearRatio];
        U1= fsolve_TE3(@fun_eep,U10,par,tol,1000);  
        U2=M1/U1;
        U3=GearRatio/M1;
        stageratio=[U1;U2;U3];
    elseif strcmp(GearConfig,'eep_2');  % speed ratio calculator when last stage ratio equals three
        Kr_1=0;
        Kr_2=1.6; %2nd stage structure weight coefficient
        M1=GearRatio/2;U10=5;
        par=[M1,Kr_1,Kr_2,Np(1:2),GearRatio];
        U1= fsolve_TE3(@fun_eep,U10,par,tol,1000);  
        U2=M1/U1;
        U3=GearRatio/M1;
        stageratio=[U1;U2;U3];
    end
    
end

