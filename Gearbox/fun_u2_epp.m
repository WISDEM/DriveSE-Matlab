
function [z]=fun_u2_epp(U,par)

M=par(1);
Kr=par(2);
Np=par(3);
M0=par(4);



z=(M^2+1)/M0/M*U^3-2*M/M0*U-2*M/M0*(M+1);



 