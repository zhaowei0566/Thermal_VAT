function dxdt=fun_expression_4_ode45(t,x,KCM)


dxdt_1 = x(2);
dxdt_2 = -KCM(2)/KCM(1)*x(2) -KCM(3)/KCM(1)*x(1) + exp(sqrt(-1)*KCM(5).*t)*KCM(4)/KCM(1);

dxdt=[dxdt_1; dxdt_2];
end