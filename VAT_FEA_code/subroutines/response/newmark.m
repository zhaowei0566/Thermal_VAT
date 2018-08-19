
function [depl, vel, accl] = newmark(sdof,R,m,k,c)
% Newmark step-by-step time integration scheme
%========================================================================%
% Written by : Siva Srinivas Kolukula (allwayzitzme(at)gmail.com)        %
% Reference : Finite Element Procedures - Bathe                          %
% Integrates the n-DOF system                                            %                                                           %
% SYNTAX :                                                               %
%           function [depl, vel, accl] = mynewmark(sdof,R,m,k,c)         %
% INPUT :                                                                %
%           sdof : System degree's of freedom                            %
%           [R]  : External applied load                                 %
%           [m]  : Assembeled Mass Matrix                                %
%           [k]  : Assembeled Stiffness MAtrix                           %
%           [c]  : Damping Matrix                                        %
%                                                                        %
% OUTPUT :                                                               %
%           depl : Displacement Response                                 %
%           vel  : Velocity                                              %
%           accl : Acceleration                                          %
%                                                                        %
% NOte: If required the time integration parameters alpha & delta can be %
%       altered                                                          %
%       Time step, initial and final time step can be changed            %
%       Initial displacements and velocities are taken as zeros          %
%========================================================================%
clc
ti = 0. ;     % initial time step
tf = 3.36 ;   % final time step
dt = 0.28 ;   % size of the time step
nt = fix((tf-ti)/dt);      % number of time steps
% initializing the displacement, velocity and acceleration matrices
depl = zeros(sdof,nt) ;
vel = zeros(sdof,nt) ;
accl = zeros(sdof,nt);
Reff = zeros(sdof,nt) ;
accl(:,1) = [0;10];
% Solve for initial accelerations
accl(:,1) = inv(m)*(R-c*vel(:,1)-k*depl(:,1));
% Parameters for Newmark time integration
alpha = 0.25 ;delta = 0.5 ;
% Calculating integration constants
a0 = 1/(alpha*dt^2) ; a1 = delta/(alpha*dt) ; a2 = 1/(alpha*dt) ;
a3 = (1/(2*alpha))-1 ; a4 = (delta/alpha)-1 ;a5 = (dt/2)*(delta/alpha-2) ;
a6 = dt*(1-delta) ; a7 = delta*dt ;
% calculating effectvie stiffness matrix
keff = k+a0*m+a1*c ;
%
% time step starts
for it = 1:nt
 Reff(:,it) = R+m*(a0*depl(:,it)+a2*vel(:,it)+a3*accl(:,it)).....
              +c*(a1*depl(:,it)+a4*vel(:,it)+a5*accl(:,it));
% solving for displacements at time (it+dt)
             depl(:,it+1)= keff\Reff(:,it);
% calculating velocities and accelerations at time (it+dt)
             accl(:,it+1) = a0*(depl(:,it+1)-depl(:,it))-a2*vel(:,it)....
             -a3*accl(:,it) ;
             vel(:,it+1) = vel(:,it)+a6*accl(:,it)+a7*accl(:,it+1);
end