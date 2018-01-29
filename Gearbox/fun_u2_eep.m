
function [z]=fun_u2_eep(U,par)

M=par(1);
Kr_1=par(2);
Kr_2=par(3);
Np_1=par(4);
Np_2=par(5);
M0=par(6);

% explict form of the dervative
% z=-M*(U-1)*(1+Kr*M*(U-1)^2)/M0/Np_2/U^2/(U/2-1)^2+2*Kr*M*(U-1)/M0/Np_2/U/(U/2-1)...
%     +(M/4/M0+Kr*M/M0/Np_2+1/M/M0)+...
%     (-M/M0/Np_2-Kr*M/M0/Np_2-M/M0)/U^2-2*M^2/M0/U^3;

% finite difference to compute the derivative, used to check the accuracy of the form above 
dU=2e-5;
v_1=M/M0/(U-dU)*(1/Np_2+1/Np_2/((U-dU)/2-1)+((U-dU)/2-1)+((U-dU)/2-1)^2+Kr_2*((U-dU)-1)^2/Np_2+Kr_2*((U-dU)-1)^2/Np_2/((U-dU)/2-1))+...
    1/M0*(1+(U-dU)/M+M/(U-dU)+M^2/(U-dU)^2);
v_2=M/M0/(U+dU)*(1/Np_2+1/Np_2/((U+dU)/2-1)+((U+dU)/2-1)+((U+dU)/2-1)^2+Kr_2*((U+dU)-1)^2/Np_2+Kr_2*((U+dU)-1)^2/Np_2/((U+dU)/2-1))+...
    1/M0*(1+(U+dU)/M+M/(U+dU)+M^2/(U+dU)^2);
z=(v_2-v_1)/2/dU;



 