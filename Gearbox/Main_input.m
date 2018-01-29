% System engineering gearbox sizing tool


% This model designs wind turbine gearboxes for minimzing the overall weight.

% Model input: rated torque (or power rating and efficiency), number of
% gearbox stages, number of planets, gearbox speed ratio, general
% configuration (eep or epp).

% Model output: gearbox weight/volume per stage, speed ratio per stage


clear all;
close all;
clc;
format short e;




dataset=[];
%input parameter for 750kW
MachineRating=750e3;
GearRatio=81.491;
GearConfig='epp';
Np=[3 1 1]; % number of planets per stage
RotorSpeed = 22;%%rpm
DrivetrainEfficiency = 0.944;
RotorTorque = (MachineRating / DrivetrainEfficiency) / (RotorSpeed * (pi / 30));

Type_ratio='optimal'; % the other choice is 'empirical'
Type_sh='normal'; % the other choice is 'short'
[Gearboxweight]=gearboxweight_estimate(GearRatio,RotorTorque,Np,GearConfig,Type_ratio,Type_sh); % get gearbox weight
dataset=[dataset;RotorTorque, Gearboxweight, MachineRating];
clear Gearboxweight RotorTorque;
disp('Gearbox weight =');
disp(dataset(:,2))


% senstivity study weigh vs eep/epp
MachineRating_s=0.5e6:0.1e6:5e6;
GearRatio=120;
Np=[4 3 1];
RotorSpeed = 15.1;
DrivetrainEfficiency = 0.95;
Type_ratio='optimal';
Type_sh='normal';
dataset_s=[];
for i=1:length(MachineRating_s)
RotorTorque_s(i) = (MachineRating_s(i)  / DrivetrainEfficiency) / (RotorSpeed * (pi / 30)); 
[Gearboxweight_epp]=gearboxweight_estimate(GearRatio,RotorTorque_s(i),Np,'epp',Type_ratio,Type_sh);
[Gearboxweight_eep]=gearboxweight_estimate(GearRatio,RotorTorque_s(i),Np,'eep',Type_ratio,Type_sh);
dataset_s=[dataset_s;RotorTorque_s(i), Gearboxweight_epp,Gearboxweight_eep, MachineRating_s(i)];
end

figure(3);plot(dataset_s(:,4),dataset_s(:,2),'-bo',dataset_s(:,4),dataset_s(:,3),'-ro');hold on;
xlabel('Rotor Torque (Nm)');ylabel('Gearbox Weight (kg)');
hleg1 = legend('Model_epp','Model_eep');set(hleg1,'Location','NorthWest');set(hleg1,'Interpreter','none')




 
