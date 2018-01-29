 function [stagemass]=stage_mass(speedratio,Np,type)
 %type=1 parallel,type=2, planetary gear;
 
 format short e;
Kr=0.4; % factor to include ring/housing/carrier weight.
if Np==3
k_gamma=1.1;
elseif Np==4
    k_gamma=1.25;
elseif Np==5
    k_gamma=1.35;
else
end
if type==1
    stagemass=1+speedratio+speedratio^2+speedratio^(-1);
elseif type==2
    sunwheelratio=0.5*speedratio-1;
    stagemass=1/Np+1/Np/sunwheelratio+sunwheelratio+sunwheelratio^2+Kr*(speedratio-1)^2/Np+Kr*(speedratio-1)^2/Np/sunwheelratio;
    stagemass=k_gamma*stagemass;
else
end
