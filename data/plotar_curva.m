clc; clear; close;
load curvaTeorica
Vmed = [IVCurve.V];   % Vetor de tensoes medidas  [V]
Imed = [IVCurve.I];
yyaxis left
plot(Vmed,Imed)
ylabel('Corrente (A)', 'FontSize', 18)
ylim([0, max(Imed)*1.1])
xlabel('Tensão (V)', 'FontSize', 18)
hold on

yyaxis right;
plot(Vmed,Imed.*Vmed)
ylim([0, max(Imed.*Vmed)*1.1])
ylabel('Potência (W)')
ylabel('Potência (W)', 'FontSize', 18)
hold on

yyaxis left;
plot(Vmed(1), Imed(1), '*k');

plot(Vmed(end), Imed(end), 'k*');
[~,id] = max(Vmed.*Imed);
plot(Vmed(id), Imed(id), '*');
vert = 0:Imed(id)/3:Imed(id);
plot(Vmed(id).*ones(size(vert)), vert, '--m');
plot(Vmed(id).*ones(size(vert)), 0, '*m');

horz = 0:Vmed(id)/3:Vmed(id);
plot(horz, Imed(id).*ones(size(horz)), '--m');
plot(0, Imed(id).*ones(size(horz)), '*m');

yyaxis right;
plot(Vmed(id), Imed(id)*Vmed(id), '*');


