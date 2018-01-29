
function [z]=fun_u1_epp(U,par)

M=par(1);
Kr=par(2);
Np=par(3);

z=(1/M^2+(Np+4*Kr)/4/Np)*U^3*(U/2-1)^2 -(Kr+Np+1)/Np*U*(U/2-1)^2 ...
    -2*M*(U/2-1)^2+2*Kr*(U-1)*U^2*(U/2-1)/Np-(Kr*(U-1)^2+1)*U*(U-1)/Np;


 