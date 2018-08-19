% function [disp,vel,accel]=harmonic_response(FEM,force,time_step,max_time)



function harmonic_response

% function second_oder_ode
%
% SOLVE  d2x/dt2+5 dx/dt - 4 x = sin(10 t)
% initial conditions: x(0) = 0, x'(0)=0

m=1;
c=5;
k=-4;
f=1;
force_omega=10;


ti=0:0.001:3;   % time scale
KCM=[m c k f force_omega];
initial_x    = 0;
initial_dxdt = 0;
xi=[initial_x initial_dxdt];
[t,x]=ode45( @fun_expression_4_ode45, ti, xi,KCM);

plot(t,x(:,1));
xlabel('t'); ylabel('x');


    function dxdt=fun_expression_4_ode45(t,x,KCM)
        
        
        dxdt_1 = x(2);
        dxdt_2 = -KCM(2)/KCM(1)*x(2) -KCM(3)/KCM(1)*x(1) + exp(sqrt(-1)*KCM(5).*t)*KCM(4)/KCM(1);
        
        dxdt=[dxdt_1; dxdt_2];
    end

end