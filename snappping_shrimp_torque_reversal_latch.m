%%%%%%%
% This is the dynamic model for the torque reversal latch modeled under the
% sets of assumptions outlined in the appendix document.
% In this model, the loading and unlatching are serial and this model goes
% through the loading, unlatching, and recoil of the spring. 
% 2022-08-01 RS

clear all; close all; clc;
%density of water
rho=1000; %kg/m^3

%design variables 
L=5.756e-3; %m - dactyl length
theta0=10*pi/180; %rad - initial dactyl opening angle (0 is fully open) %defined slightly different than previous model for convenience
D=1.025e-3; %m - estimate of diameter cylinder for mass estimate
m=(pi/4*L*D^2)*rho; %kg - mass of dactyl (volume - m^3 * 1000 kg/m^3)
% this is setting the mass of the dactyl similar to that of an equivalent
% volume of water
Icm=1/12*m*L^2; %kg-m^2 mass moment of inertia at the center of a bar
%calculated from design
I=Icm+m*((L/2)^2+(D/2)^2); %kg-m^2 mass moment of inertia at the pivot 

%this is used to calculate the moment of inertia about the pivot point
av=.264e-3*1; %m - vertical offset from pivot to tendon
ah=D; %m - horizontal offset from pivot to tendon

% calculating the trigger angle in degrees and radians
GeoAngle=atand(av/ah); % deg
theta_tr=atan(av/ah); % rad

%spring constant
kt=5000; %N/m - elastic element spring constant

% drag coefficients
Cd=0.5;



% set up time step and simulation time
dt=1e-6;
time=[0:dt:1]';

%input energy to system - this will be setting how much energy is desired
%overall, but the dynamics will take care of if that energy is achieved or
%not
E0=6e-3; %J - desired elastic energy to system
t0=100e-3*4; %s - desired time for elastic energy to be loaded 
% (energy loading duration)
%the ratio of these two is the energy loading rate


Evec(1,1)=0; % initial point to vector the for loop will step through this 
% loading profile in time. We could also change how this is done here in
% this step, right now it is just a linear fit that maxes out
for indexenergy=2:length(time)
    Evec(indexenergy,1)=min((Evec(indexenergy-1,1)+dt*E0/t0),E0); 
    %building energy vector in relation to simulation time
end

% calculating the deformation of the spring
xt=sqrt(2*Evec./kt);


% loading of spring
% assumption for now, loading is separted from the tendon moving down,
% though this can adjusted such that loading can continue while the spring
% is shifted down

dtheta(1,1)=0;
theta(1,1)=theta0;

% loading sequence of mechanism
for indexl=1:1:(t0/dt)
    xnew(indexl,1)=sqrt(xt(indexl,1).^2-I/kt*dtheta(indexl,1).^2); % this
    % allows for the elastic recoil of the spring without monitoring both
    % ends of the spring. In this we are subtracting the kinetic energy of
    % the dactyl movement from the spring potential energy
    Ft(indexl,1)=kt*xnew(indexl,1); %calculates the new elastic force
    theta_ul(indexl,1)=0; %rad - spring actuation angle
    
    %checking where the dactyl motion and unlatching angle is relative to
    %the critical angle. The dynamics and moment arm changes after
    %over-centering
    %alpha is placeholder variable
    if (theta(indexl,1)+theta_ul(indexl)) < theta_tr
        alpha(indexl,1)=theta_tr-(theta(indexl,1)+theta_ul(indexl,1));
        Torque_dist(indexl,1)=-1*sqrt(ah^2+av^2)*tan(alpha(indexl,1));
    else
        alpha(indexl,1)=(theta(indexl,1)+theta_ul(indexl,1))-theta_ul(indexl,1);
        Torque_dist(indexl,1)=1*sqrt(ah^2+av^2)*sin(alpha(indexl,1));
    end
    
    ddtheta(indexl,1)=(Ft(indexl,1)*Torque_dist(indexl,1)-...
        Cd*rho*D*L^4/8*dtheta(indexl,1)^2)/I;

    %in this statement, if the torque is acting in the "wrong" direction
    %this statement will similar a mechanical block by not allowing motion
    if ddtheta(indexl,1)<0
        dtheta(indexl+1,1)=0*(dtheta(indexl,1)+ddtheta(indexl,1)*dt);
        theta(indexl+1,1)=1*(theta(indexl,1)+dtheta(indexl,1)*dt);
    else
        dtheta(indexl+1,1)=dtheta(indexl,1)+ddtheta(indexl,1)*dt;
        theta(indexl+1,1)=theta(indexl,1)+dtheta(indexl,1)*dt;
    end
    
    %code ends when claw closes
    if abs(theta(indexl,1))>pi/2
        break
    end
end

%unlatching rate
dtheta_ul_dt=1 * (45*pi/180)/.1; %rad/s - rotation rate of tendon to 
% over centering

%variable to catch when over-centering happens
unlatchingindex=0;

% dynamics and code in this for loop are the same as above
for indexul=indexl+1:length(time)
    xnew(indexul,1)=sqrt(xt(indexul,1).^2-I/kt*dtheta(indexul,1).^2); %calculates elastic recoil
    Ft(indexul,1)=kt*xnew(indexul,1);
    theta_ul(indexul,1)=theta_ul(indexul-1,1)+dtheta_ul_dt*dt;
    
    if (theta(indexul,1)+theta_ul(indexul)) < theta_tr
        alpha(indexul,1)=theta_tr-(theta(indexul,1)+theta_ul(indexul));
        Torque_dist(indexul,1)=-1*sqrt(ah^2+av^2)*tan(alpha(indexul,1)); %placeholder for distance
    else
        alpha(indexul,1)=(theta(indexul,1)+theta_ul(indexul))-theta_ul(indexl,1);
        Torque_dist(indexul,1)=1*sqrt(ah^2+av^2)*sin(alpha(indexul,1));
    end
    
    ddtheta(indexul,1)=(Ft(indexul,1)*Torque_dist(indexul,1)-...
        Cd*rho*D*L^4/8*dtheta(indexul,1)^2)/I;
    
    if ddtheta(indexul,1)<0
        dtheta(indexul+1,1)=0*(dtheta(indexul,1)+ddtheta(indexul,1)*dt);
        theta(indexul+1,1)=1*(theta(indexul,1)+dtheta(indexul,1)*dt);
    else
        dtheta(indexul+1,1)=dtheta(indexul,1)+ddtheta(indexul,1)*dt;
        theta(indexul+1,1)=theta(indexul,1)+dtheta(indexul,1)*dt;
        if unlatchingindex==0
            unlatchingindex=indexul;
        end
    end
    

        if abs(theta(indexul,1))>=(pi/2)
           break
        end
end

strikeduration=time(indexul,1)-time(unlatchingindex,1);
strikeclosure=theta(indexul,1)-theta(indexl,1);
