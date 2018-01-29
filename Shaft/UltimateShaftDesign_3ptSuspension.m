%NREL NWTC
%10/21/2013
%

%This script performs preliminary main shaft design for tip speed study

%The script is coded based on a analytical model that outputs drivetrain loads and bending moments.
%The model assumes the gearbox as a solid piece.

%The program input is the aerodynamic forces at hub center, both modeling
%and experimental results can be used.

%The program output is the maximum shaft length of LSS

%All the results shown below are reaction forces in metric units
format short e;
clear all;
close all;
clc;

% load load_baseline_100ms_flexible;
load load_baseline_750kW;

m_hub=110000.234*0;
L_overhung=2.3758;
MachineRating=0.75e6;
RotorSpeed = 22.0;
DrivetrainEfficiency = 0.944;
RotorTorque = (MachineRating  / DrivetrainEfficiency) / (RotorSpeed * (pi / 30))/1e3;

%Input Parameters
g=9.81;         % Gravity G with consideratin of 6 degree tilting angle
gamma=5;        % Tilting angle
PSF=1;

tol=1E-4;




index=[1:12];
for iii=1:1%length(loads_ms(:,1))  % deterministic cases 9, 10, 11, 12, 21
    close all
    clear x_shaft
    clear My_ms
    clear Mz_ms
    clear theta_y
    clear d_y
    counter=0;
    check_limit=1;
    L_ms_new=0;
    L_ms_0=0.5;     % Main shaft length
    dL=0.05;
    L_ms=L_ms_0;
    D_max=1;
    D_min=0.2;
    while abs(check_limit)>tol && counter<50
        counter=counter+1;
        if L_ms_new>0
            
            L_ms=L_ms_new;
        else
            L_ms=L_ms_0;
        end
        %Distances
        L_rb=1.912;%L_overhung;%1.912;     % Distance between hub center and MB1;
        L_bg=6.11;      % Distance between hub center and gearbox yokes
        
        L_as=L_ms/2;   % Distance between main bearing and shaft center
        L_gb=0;         % Gearbox center w.r.t. gearbox trunnions in the x direction
        H_gb=1;         % Gearbox center w.r.t. gearbox trunnions in the z direction
        L_gp=0.825;     % Distance between gearbox coupling and gearbox trunnions
        L_cu=L_ms+0.5;
        L_cd=L_cu+0.5;
        %Weight
        W_r=m_hub*g;   % Weight of hub and blades
        W_ms=pi/3*(D_max^2+D_min^2+D_max*D_min)*L_ms*7800*g/4;  % Weight of main shaft
        
        W_a=1000*g;  % Addtional mass at the main shaft
        W_gb=5.3692e4*g/5;   % Gearbox weight
        W_c=1e3*g;  % Carrier Weight
        %Hub forces
        F_r_x=loads_ms(iii,1)*1e3/PSF;  % External F_x from the blades
        F_r_y=loads_ms(iii,2)*1e3/PSF;  % External F_y from the blades
        F_r_z=loads_ms(iii,3)*1e3/PSF;  % External F_z from the blades
        M_r_x=loads_ms(iii,4)*1e3/PSF;  % External M_x from the blades
        T=loads_ms(iii,4)/PSF;%4.9181E+03;%RotorTorque*1.25/20;% kNm
        M_r_y=loads_ms(iii,5)*1e3/PSF;  % Nm External M_y from the blades
        M_r_z=loads_ms(iii,6)*1e3/PSF;  % External M_z from the blades
        % material properties
        E=2.1e11;
        density=7850;
        n_safety=2.5;%1.93;
        Sy=66000; %psi
        % OD at main bearing
        u_knm_inlb=8850.745454036;
        u_in_m=0.0254000508001;
        %bearing deflection check
        MB_limit=0.026;
        CB_limit=4/60/180*pi;
        TRB_limit=3/60/180*pi;
        n_safety_brg=1.0;
        
        %define main shaft
        len_pts=101;
        x_ms=linspace(L_rb,L_ms+L_rb,len_pts);
        x_rb=linspace(0,L_rb,len_pts);
        y_gp=linspace(0,L_gp,len_pts);
        len=1:length(M_r_y);
        
        for i=1:len
            
            f_mb_x(i)=-F_r_x(i)-W_r*sind(gamma);
            f_mb_y(i)=M_r_z(i)/L_bg-F_r_y(i)*(L_bg+L_rb)/L_bg;
            f_mb_z(i)=(-M_r_y(i)+W_r*(cosd(gamma)*(L_rb+L_bg)+sind(gamma)*H_gb)+...
                W_ms*(L_bg-L_as)*cosd(gamma)+W_a*cosd(gamma)*(L_bg-L_ms)+...
                -W_gb*cosd(gamma)*L_gb-F_r_z(i)*cosd(gamma)*(L_bg+L_rb))...
                /L_bg;
            
            
            f_gb_x(i)=-(W_ms+W_a+W_gb)*sind(gamma);
            f_gb_y(i)=-f_mb_y(i)-F_r_y(i);
            f_gb_z(i)=-f_mb_z(i)+(W_a+W_r+W_ms+W_gb)*cosd(gamma)-F_r_z(i);
            
            f_cu_z(i)=(W_ms*cosd(gamma)+W_a*cosd(gamma)+W_gb*cosd(gamma))-f_mb_z(i)-F_r_z(i)-...
                (-M_r_y(i)-F_r_z(i)*cosd(gamma)*L_rb+W_ms*(L_bg-L_as)*cosd(gamma)-W_c*cosd(gamma)*L_gb)/(1-L_cu/L_cd);
            f_cd_z(i)=(W_ms*cosd(gamma)+W_a*cosd(gamma)+W_gb*cosd(gamma))-f_mb_z(i)-F_r_z(i)-f_cu_z(i);
            
            for k=1:length(x_rb)
                My_ms(k)=-M_r_y(i)+W_r*cosd(gamma)*x_rb(k)+0.5*W_ms/L_ms*x_rb(k)^2-F_r_z(i)*x_rb(k);
                Mz_ms(k)=-M_r_z(i)-F_r_y(i)*x_rb(k);
            end
            
            for j=1:length(x_ms)
                My_ms(j+len_pts)=-F_r_z(i)*x_ms(j)-M_r_y(i)+W_r*cosd(gamma)*x_ms(j)-f_mb_z(i)*(x_ms(j)-L_rb)+0.5*W_ms/L_ms*x_ms(j)^2;
                Mz_ms(j+len_pts)=-M_r_z(i)-f_mb_y(i)*(x_ms(j)-L_rb)-F_r_y(i)*x_ms(j);
            end
            
        end
        
        x_shaft=[x_rb,x_ms];
        figure(1);subplot(2,1,1);plot(x_shaft,My_ms/1e3,'linewidth',2);hold on;plot(x_shaft,Mz_ms/1e3,'r','linewidth',2);
        hold on;plot(x_shaft,(My_ms.^2+Mz_ms.^2).^0.5/1e3,'k--','linewidth',2);
        hold on;plot(x_rb(end),(My_ms(length(x_rb)))/1e3,'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','r',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',8)
        xlabel('Main Shaft Length, m','fontsize',12);ylabel('Shaft Bending, kNm','fontsize',12);
        text(x_rb(end),1.1*(My_ms(length(x_rb)))/1e3,['Main Bearing'],'FontSize',10);hold on;
        hleg = legend('My','Mz','Norm',...
            'Location','SouthEast');
        axis fill;
        
        [MM_max,Index]=max((My_ms.^2+Mz_ms.^2).^0.5/1e3);
        disp('Amplitude = kNm');disp(MM_max);disp('Location = m');disp(x_shaft(Index));
        MM_min=((My_ms(end).^2+Mz_ms(end).^2).^0.5/1e3);
        disp('Amplitude = kNm');disp(MM_min);disp('Location = m');disp(x_shaft(end));
        
        %%%%%%%%%%%%%%%%
        %design shaft OD
        
        MM=MM_max;
        D_max=(16*n_safety/pi/Sy*(4*(MM*u_knm_inlb)^2+3*(T*u_knm_inlb)^2)^0.5)^(1/3)*u_in_m;
        % OD at end
        MM=MM_min;
        D_min=(16*n_safety/pi/Sy*(4*(MM*u_knm_inlb)^2+3*(T*u_knm_inlb)^2)^0.5)^(1/3)*u_in_m;
        
        
        % Estimate ID
        D_in=0.10*D_max;
        D_max=(D_in^4+D_max^4)^0.25;
        D_min=(D_in^4+D_min^4)^0.25;
        disp('shaft OD: max = m')
        disp(D_max)
        disp('shaft OD: min = m')
        disp(D_min)
        
        % Plot MS
        x_sh_plot=[0,0,x_rb(end),x_ms(end),x_ms(end),x_rb(end),0];
        y_sh_plot=[-D_max/2,D_max/2,D_max/2,D_min/2,-D_min/2,-D_max/2,-D_max/2];
        subplot(2,1,2);plot(x_sh_plot,y_sh_plot,'k','linewidth',2); axis fill;hold on;
        xlabel('Main Shaft Length, m','fontsize',12);ylabel('Main Shaft OD, m','fontsize',12);
        x_sh_plot_in=[0,x_ms(end)];y_sh_plot_in_up=[D_in/2,D_in/2];y_sh_plot_in_dn=[-D_in/2,-D_in/2];
        plot(x_sh_plot_in,y_sh_plot_in_up,'k--','linewidth',2);plot(x_sh_plot_in,y_sh_plot_in_dn,'k--','linewidth',2)
        hold on;plot(x_rb(end),0,'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','r',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',8);
        
        
        W_ms_new=(density*pi/12*L_ms*(D_max^2+D_min^2+D_max*D_min)...
            -density*pi/4*D_in^2*L_ms...
            +density*pi/4*D_max^2*L_rb)*g;
        
        
        fx=vectorize(inline('-F_r_z*z^3/6+W_r*cosd(gamma)*z^3/6-M_r_y*z^2/2-f_mb_z*(z-L_rb)^3/6+W_ms/(L_ms+L_rb)/24*z^4','F_r_z','W_r','gamma','M_r_y','f_mb_z','L_rb','W_ms','L_ms','z'));
        D1=fx(F_r_z,W_r,gamma,M_r_y,f_mb_z,L_rb,W_ms_new,L_ms,L_rb+L_ms);
        D2=fx(F_r_z,W_r,gamma,M_r_y,f_mb_z,L_rb,W_ms_new,L_ms,L_rb);
        C1=-(D1-D2)/L_ms;
        C2=-D2-C1*(L_rb);
        
        I_2=pi/64*(D_max^4-D_in^4);
        gx=vectorize(inline('-F_r_z*z^2/2+W_r*cosd(gamma)*z^2/2-M_r_y*z-f_mb_z*(z-L_rb)^2/2+W_ms/(L_ms+L_rb)/6*z^3+C1','F_r_z','W_r','gamma','M_r_y','f_mb_z','L_rb','W_ms','L_ms','C1','z'));
        
        for kk=1:len_pts
            theta_y(kk)=gx(F_r_z,W_r,gamma,M_r_y,f_mb_z,L_rb,W_ms_new,L_ms,C1,x_ms(kk))/E/I_2;
            d_y(kk)=(fx(F_r_z,W_r,gamma,M_r_y,f_mb_z,L_rb,W_ms_new,L_ms,x_ms(kk))+C1*x_ms(kk)+C2)/E/I_2;
        end
        
        figure(5);plot(x_ms,theta_y,'k');
        hold on;plot(x_rb(end),(theta_y(1)),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','r',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',8)
        xlabel('Main Shaft Length, m','fontsize',12);ylabel('Shaft Angular Deflection, rad','fontsize',12);
        text(x_rb(end),1.1*(theta_y(1)),['Main Bearing'],'FontSize',10);hold on;
        % plot(x_ms,ones(length(x_ms))*MB_limit,'r+');text(x_rb(end),1.1*(MB_limit),['Main Bearing limit'],'FontSize',10);hold on;
        plot(x_ms,ones(length(x_ms))*TRB_limit,'r+');text(x_ms(end)*0.9,1.0*(TRB_limit),['Carrier Bearing limit'],'FontSize',10);hold on;
        hold on;plot(x_ms(end),(theta_y(end)),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','r',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',8)
        text(x_ms(end),1.1*(theta_y(end)),['Carrier Bearing'],'FontSize',10);hold on;axis fill;
        plot(x_ms,-ones(length(x_ms))*TRB_limit,'r+');text(x_ms(end)*0.9,1.0*(TRB_limit),['Carrier Bearing limit'],'FontSize',10);hold on;
        hold on;plot(x_ms(end),(theta_y(end)),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','r',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',8)
        text(x_ms(end),1.1*(theta_y(end)),['Carrier Bearing'],'FontSize',10);hold on;axis fill;
        % L_ms_new=max(x_ms(find(abs(theta_y-TRB_limit/n_safety_brg)<1e-5))-L_rb);
        
        
        check_limit=abs(abs(theta_y(end))-TRB_limit/n_safety_brg);
        if check_limit<0
            L_ms_new=L_ms+dL;
        else
            L_ms_new=L_ms+dL;
        end
        
        disp('new shaft length = m');
        disp(L_ms_new);
        % figure;plot(x_ms,d_y)
        pause(0.05)
    end
end

%%
 %%% plot
    
     x_sh_plot=[x_rb(end)-0.2,x_rb(end)-0.2,x_rb(end)+0.2,x_mb(end)-0.25,x_mb(end)+0.25,x_mb(end)+0.25,x_mb(end)-0.25,x_rb(end)+0.2,x_rb(end)-0.2];
            y_sh_plot=[-D_max/2,D_max/2,D_max/2,D_medium/2,D_medium/2,-D_medium/2,-D_medium/2,-D_max/2,-D_max/2];
            subplot(2,1,2);plot(x_sh_plot,y_sh_plot,'k','linewidth',2); axis fill;hold on;
            xlabel('Main Shaft Length, m','fontsize',12);ylabel('Main Shaft OD, m','fontsize',12);
            x_sh_plot_in=[x_rb(end)-0.2,x_mb(end)+0.25];y_sh_plot_in_up=[D_in/2,D_in/2];y_sh_plot_in_dn=[-D_in/2,-D_in/2];
            plot(x_sh_plot_in,y_sh_plot_in_up,'k--','linewidth',2);plot(x_sh_plot_in,y_sh_plot_in_dn,'k--','linewidth',2)
            x_sh_plot_vertical_2=[x_mb(end)-0.25,x_mb(end)-0.25];y_sh_plot_in_vertical_2=[-D_medium/2,D_medium/2];
            plot(x_sh_plot_vertical_2,y_sh_plot_in_vertical_2,'k--','linewidth',2);
            x_sh_plot_vertical_1=[x_rb(end)+0.2,x_rb(end)+0.2];y_sh_plot_in_vertical_1=[-D_max/2,D_max/2];
            plot(x_sh_plot_vertical_1,y_sh_plot_in_vertical_1,'k--','linewidth',2)
            
            hold on;plot(x_rb(end),0,'-mo',...
                'LineWidth',2,...
                'MarkerEdgeColor','r',...
                'MarkerFaceColor',[.49 1 .63],...
                'MarkerSize',8);
            hold on;plot(x_mb(end),0,'-mo',...
                'LineWidth',2,...
                'MarkerEdgeColor','r',...
                'MarkerFaceColor',[.49 1 .63],...
                'MarkerSize',8);
            density=7800;
            W_ms_new=(density*pi/12*L_ms*(D_max^2+D_min^2+D_max*D_min)...
                -density*pi/4*D_in^2*L_ms...
                +density*pi/4*D_max^2*L_rb)*g;
            
%             D_max=1.0;D_medium=0.63; L_mb=2.80-0.175-0.2;%4.434-L_rb+0.2;
%             total_weight=pi/3*(D_max^2+D_medium^2+D_max*D_medium)*(L_mb)*7800/4+...
%                          pi/4*(D_max^2-D_in^2)*7800*0.35+pi/4*(D_medium^2-D_in^2)*7800*0.4-...
%                          pi/4*(D_in^2)*7800*(L_mb+0.75/2);
%                       
            D_max=0.35;D_medium=0.34; L_mb=1.4;
            total_weight=pi/3*(D_max^2+D_medium^2+D_max*D_medium)*(L_mb)*7800/4+...
                         pi/4*(D_max^2-D_in^2)*7800*0.18+pi/4*(D_medium^2-D_in^2)*7800*0.18-...
                         pi/4*(D_in^2)*7800*(L_mb+0.18*2);

    
    