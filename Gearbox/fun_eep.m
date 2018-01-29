
function [z]=fun_eep(U,par)

M=par(1);
Kr_1=par(2);
Kr_2=par(3);
Np_1=par(4);
Np_2=par(5);
M0=par(6);
% explict form of the dervative

% z=-(U-1)*(1+Kr*(U-1)^2)/Np_1/U^2/(U/2-1)^2+2*Kr*(U-1)/Np_1/U/(U/2-1)-(1+Kr)/Np_1/U^2+0.25+Kr/Np_1...
%      +0.5/Np_2/(M/2-U)^2-2*Kr*(M-U)/Np_2/U^2/(M/2-U)+Kr*(M-U)^2/2/Np_2/U^2/(M/2-U)^2 ...
%      +(0.5+2*Kr/Np_2)*(1/U^2-M/U^3);

% finite difference to compute the derivative, used to check the accuracy
% of the form above

dU=2e-5;
v_1=1/(U-dU)*(1/Np_1+1/Np_1/((U-dU)/2-1)+((U-dU)/2-1)+((U-dU)/2-1)^2+Kr_1*((U-dU)-1)^2/Np_1+Kr_1*((U-dU)-1)^2/Np_1/((U-dU)/2-1))+...
1/M*(1/Np_2+1/Np_2/(M/2/(U-dU)-1)+(M/2/(U-dU)-1)+(M/2/(U-dU)-1)^2+Kr_2*(M/(U-dU)-1)^2/Np_2+Kr_2*(M/(U-dU)-1)^2/Np_2/(M/2/(U-dU)-1));
v_2=1/(U+dU)*(1/Np_1+1/Np_1/((U+dU)/2-1)+((U+dU)/2-1)+((U+dU)/2-1)^2+Kr_1*((U+dU)-1)^2/Np_1+Kr_1*((U+dU)-1)^2/Np_1/((U+dU)/2-1))+...
1/M*(1/Np_2+1/Np_2/(M/2/(U+dU)-1)+(M/2/(U+dU)-1)+(M/2/(U+dU)-1)^2+Kr_2*(M/(U+dU)-1)^2/Np_2+Kr_2*(M/(U+dU)-1)^2/Np_2/(M/2/(U+dU)-1));
z=(v_2-v_1)/2/dU;