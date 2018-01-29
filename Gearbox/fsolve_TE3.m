function [x1]=fsolve_TE3(fhandle,x0,parameters,tol,max_fsolve_iter)
global fsolve_iter

x0=x0.';
fsolve_iter=0;

error=10;

n=length(x0);
while (fsolve_iter<max_fsolve_iter) && (error>tol);
    fsolve_iter=fsolve_iter+1;
    for i=1:n
        e(i,:)=[x0(i)*0.95 x0(i)*1.05];
        h(i)=e(i,2)-e(i,1);
        xplus=x0;
        xminus=x0;
        
        xplus(i)=e(i,2);
        xminus(i)=e(i,1);
        
        Jacobian(:,i)=(fhandle(xplus,parameters)-fhandle(xminus,parameters))/h(i);
    end
    f=fhandle(x0,parameters).';


    error=sqrt(sum(f.^2));
    x1=x0-Jacobian\f;
    x0=x1;
end
if fsolve_iter==max_fsolve_iter 
    disp('Max fsolve_iterations reached');
    x1=NaN+zeros(size(x1));
end

if abs(sum(imag(x0)))>0
    disp('Imaginary');
    x1=NaN+zeros(size(x1));
end