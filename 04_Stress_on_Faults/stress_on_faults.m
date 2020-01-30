% This script calculates the shear and normal stress on faults of any
% orientation given a stress state and pore pressure. It also calculates
% the likelyhood of slip using the Coulomb failure criterion. 

function [FaultNo_Strike_Dip_ShearStress_EffNormalStress]=stress_on_faults(S1,S2,S3,Pp,mu,a,b,c,faults) 

% S1, S2, S3 = Magnitude of principal stresses (MPa)
% Pp = Pore pressure (MPa)
% mu = Coefficient of sliding friction
% a = trend of S1, except when S1 is vertical a = trend of SHmax minus 90 degrees
% b = -plunge of S1 (plunge is angle from horizontal)
% c = rake of S2, 0 if S1 or S3 is vertical, 90 if S2 is vertical
% faults = 2 column matrix with strike(0-360), dip(0-90) 
%          (use right-hand rule)), e.g faults=[208,60;36,62;...]
%               
% All angle information should be input as degrees, stresses as MPa (output
% will be the same

format short

% Converting fault and stress orientation info from degrees to radians 
str=deg2rad(faults(:,1));
dip=deg2rad(faults(:,2));
a=deg2rad(a);
b=deg2rad(b);
c=deg2rad(c);

% Define principal stress tensor
S=[S1 0 0
    0 S2 0
    0 0 S3];

% Transformation from principal stress to geographic coordinates
R1=[cos(a)*cos(b) sin(a)*cos(b) -sin(b);
    cos(a)*sin(b)*sin(c)-sin(a)*cos(c) sin(a)*sin(b)*sin(c)+cos(a)*cos(c) cos(b)*sin(c);
    cos(a)*sin(b)*cos(c)+sin(a)*sin(c) sin(a)*sin(b)*cos(c)-cos(a)*sin(c) cos(b)*cos(c)];

Sg=R1'*S*R1;

% For loop makes calculations for each fault one at a time
fault_no=[1:1:length(faults(:,1))]';
for k=1:length(fault_no)
    % Transformation from geographic to fault coordinate system
    R2=[cos(str(k)) sin(str(k)) 0;
         sin(str(k))*cos(dip(k)) -cos(str(k))*cos(dip(k)) -sin(dip(k));
         -sin(str(k))*sin(dip(k)) cos(str(k))*sin(dip(k)) -cos(dip(k))];

    Sf=R2*Sg*R2';
    
    % Solve for normal stress resolved on the fault surface
    Sn(k)=Sf(3,3);

    % Solve for rake of the slip vector
    if (Sf(3,2)>0)
        rake(k)=(atan(Sf(3,2)/(Sf(3,1))));
    elseif (Sf(3,2)<0)&(Sf(3,1)>0)
        rake(k)=pi-atan(Sf(3,2)/(-Sf(3,1)));
    else
        rake(k)=atan((-Sf(3,2))/(-Sf(3,1)))-pi;
    end
    
    % Solve for shear stress resolved on the fault plane
    R3=[cos(rake(k)) sin(rake(k)) 0;
        -sin(rake(k)) cos(rake(k)) 0;
        0 0 1];

    Sr=R3*Sf*R3';

    tau(k,1)=Sr(3,1);
    
    % Solve for Coulomb failure function (CFF) (proxy for likelihood of slip)
    Sn_eff(k,1)=Sn(k)-Pp;
    CFF(k,1)=abs(tau(k)) - mu*Sn_eff(k);
    rake(k)=rad2deg(rake(k));
end

%Output data changing signs to make output easier to understand
rake=(sign(tau).*rake');
tau=abs(tau);

%Plot CFF values on a Mohr diagram

close all

figure

beta = linspace(0,2*pi,360);

sigman = linspace(0,30);

friction = mu.*sigman;

sigma1 = S1-Pp;
sigma2 = S2-Pp;
sigma3 = S3-Pp;

x1 = ((sigma1-sigma3)/2).*cos(beta) + (sigma1+sigma3)/2;
y1 = ((sigma1-sigma3)/2).*sin(beta);

x2 = ((sigma1-sigma2)/2).*cos(beta) + (sigma1+sigma2)/2;
y2 = ((sigma1-sigma2)/2).*sin(beta);

x3 = ((sigma2-sigma3)/2).*cos(beta) + (sigma2+sigma3)/2;
y3 = ((sigma2-sigma3)/2).*sin(beta);

plot(x1,y1,'-k')

hold on

plot(x2,y2,'-k')

plot(x3,y3,'-k')

plot(sigman,friction,'-k')

scatter(Sn_eff,tau,20,CFF,'filled')

axis('equal')

axis([0 30 0 15])

xlabel('Normal stress (MPa)')
ylabel('Shear stress (MPa)')

c=colorbar;
ylabel(c,'Coulomb Failure Function (CFF)')

FaultNo_Strike_Dip_ShearStress_EffNormalStress=[fault_no,faults,tau,Sn_eff,CFF];



