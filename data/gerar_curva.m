clc; clear; close;
k = 1.3806503e-23; q = 1.60217646e-19;
Tc = 33;
Ns = 1;
T = Tc + 273.15; % temperatura [K]
Vt = Ns*k*T/q;

Iph = 0.5;
I0 = 0.5*10^-6;
n = 1;
Rs = 0.01;
Rp = 50;
x = [0.5, 0.5*10^-6, 1, 0.01, 50];

Rsh = Rp;
modelo = @(V,I) (I - Iph + I0*(exp((V+I*Rs)/(n*Vt))-1) + (V+I*Rs)/Rp);
valores = fimplicit(modelo,[0 0.38 0 Iph*1.1]);

IVCurve.V = [valores.XData].';
V = IVCurve.V;
theta = Rs*I0*Rsh * exp(Rsh*(Rs*(Iph+I0) + V)/(n*Vt*(Rs + Rsh))) / (n*Vt*(Rs + Rsh));
IVCurve.I = Rsh*(Iph + I0)/(Rs + Rsh) - V/(Rs + Rsh) - (n*Vt/Rs)*lambertw(theta);
IVCurve.Ns = Ns;
IVCurve.T = Tc;
figure
plot(IVCurve.V, IVCurve.I);
save('curvaTeorica', 'IVCurve');

RMSE2 = sqrt(fobj2(x,IVCurve.V,IVCurve.I, Vt, 1))
RMSE1 = sqrt(fobj(x,IVCurve.V,IVCurve.I, Vt, 1))
%%
function MSE = fobj2(x, Vmed, Imed, Vth, pop_size)
V = Vmed;
MSE = zeros(0,pop_size);
for i = 1:pop_size
    Iph = x(i,1);
    I0 = x(i,2);
    n = x(i,3);
    Rs = x(i,4);
    Rsh = x(i,5);
    theta = Rs*I0*Rsh * exp(Rsh*(Rs*(Iph+I0) + V)/(n*Vth*(Rs + Rsh))) / (n*Vth*(Rs + Rsh));
    Imod = Rsh*(Iph + I0)/(Rs + Rsh) - V/(Rs + Rsh) - (n*Vth/Rs)*lambertw(theta);
    MSE(i) = sum((Imod - Imed).^2)/length(Imed);
end
end

%% funcao objetivo, modelo de um diodo
function [f] = fobj(x, Vmed, Imed, Vt, pop_size)
% sqrt foi removido para melhorar performace
MSE = zeros(1,pop_size); % pre-alocacao de memoria
l = length(Imed);
for i = 1:pop_size
    Iph = x(i,1);
    I0 = x(i,2);
    n = x(i,3);
    Rs = x(i,4);
    Rp = x(i,5);
    MSE(i) = sum((Imed - Iph + I0*(exp((Vmed+Imed*Rs)/(n*Vt))-1) + (Vmed+Imed*Rs)/Rp).^2)/l;
end
f = MSE;
end
