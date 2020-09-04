load rtcFrance
hold off
T = IVCurve.T;
Ns = 1;

k = 1.3806503e-23;
q = 1.60217646e-19;
T = T + 273.15;
Vt = Ns * k * T / q;

Iph = 0.76078;
I0 = 0.32302*10^-6;
Rs = 0.03638;
Rp = 53.71636;
n = 1.48118;

% Iph = 0.5;
% I0 = 0.5*10^-6;
% n = 1;
% Rs = 0.01;
% Rp = 50;
% V = 0:0.01:0.37;

V= IVCurve.V;
I = IVCurve.I;
I_autores = Iph - I0.*(exp((V + I.*Rs)./(n*Vt))-1) - (V + I.*Rs)./Rp;

Inovo = getI(Rs, Rp, n, Iph, I0, T, Ns, V, I,Vt,20);
plot(IVCurve.V, I_autores,'.m')
hold on
plot(IVCurve.V, Inovo,'.b')

function [I, V] = getI(Rs, Rp, n, Iph, I0, T, Ns,V, Imed, Vt, iter)
% Estimativa inicial
% I = Iph - I0.*(exp((V)./(n*Vt))-1);
I = Imed;
    plot(V, I)
    drawnow
        hold on
% melhoramento iterativo
for i = 1:iter
    In = Iph - I0.*(exp((V + I.*Rs)./(n*Vt))-1) - (V + I.*Rs)./Rp;
    plot(V, I,'.')
    err(i) = max(abs(In-I))
    I=In;
    drawnow
end
end