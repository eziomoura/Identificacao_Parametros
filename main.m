clc; clear; close all;
% OBS2: CIABC foi perda de tempo. Ver notas no arquivo
% OBS: RMSE deve ser retornado como vetor coluna
%%
CODE_FUN_OBJ = 1; % ver arquivo makeFunObj, para lista de codigos
selAlgo = {'TLABC', 'TLBO'}; % Vetor com os algoritmos que deseja avaliar
listAlgo = {'BFS','ABC','DE','EJADE','IJAYA','ITLBO','JADE','PGJAYA','PSO','TLBO'}; % (nao atualizada) Lista de todos algoritmos disponíveis
RUNS = 30; % quantidade de execuções distintas
POP = 50; % tamanho da população (>5)
MAX_FES = 50000; % numero maximo de avalicoes da funcao objetivo
graphic = false; % deseja plotar curvas IV?

%% Dados de entrada
addpath('.\Outros Algoritmos')
load 'data/photowatt.mat'
plotCurvesEveryRun = false;
Vmed = [IVCurve.V];   % Vetor de tensoes medidas  [V]
Imed = [IVCurve.I];   % Vetor de correntes medidas [A]
Tc =   [IVCurve.T];   % Temperatura [ºC]
Ns =   [IVCurve.Ns];  % Numero de celulas em serie
%% Criterios de ajuste
%          var = [Iph, I0, n, Rs, Rp];
if Ns > 1
    limite_inf = [0, 0, 1, 0, 0];            % limite inferior do módulo
    limite_sup = [2, 50e-6, 50/Ns, 2, 2000]; % limite superior
else
    limite_inf = [0, 0, 1, 0, 0];            % limite inferior da célula
    limite_sup = [1, 1e-6, 2, 0.5, 100];     % limite superior
end
%% Define a função objetivo
addpath('.\Funções Objetivo')
fun = makeFobj(Vmed, Imed, Ns, Tc, CODE_FUN_OBJ);
fobj = @fun.Objective;

%% pre alocação
%converg_curve = zeros(RUNS, maxIter);
elapsedTime = zeros(1,RUNS);
Iph = zeros(RUNS,1); I0 = zeros(RUNS,1); n = zeros(RUNS,1);
Rs = zeros(RUNS,1);  Rp = zeros(RUNS,1); RMSE = zeros(RUNS,1);
if strcmp(selAlgo{1}, 'all')
    selAlgo = listAlgo;
end
for i = 1: length(selAlgo)
    metaheuristic = getAlgo(selAlgo{i}); 
    fprintf('\nTestando %s \n', selAlgo{i});
    
    for run = 1:RUNS
        fprintf('\nExecução %d \n', run)
        tic
        [x, RMSE(run), converg(run).RMSE, converg(run).fes] =  metaheuristic(fobj, limite_inf, limite_sup, POP, MAX_FES, true);
        
        Iph(run) = x(1);
        I0(run) = x(2);
        n(run) = x(3);
        Rs(run) = x(4);
        Rp(run) = x(5);
        
        elapsedTime(run) = toc;
        fprintf(' Iph(A) = %f', Iph(run));
        fprintf('\n I0(uA) = %f', I0(run)*10^6);
        fprintf('\n n = %f', n(run));
        fprintf('\n Rs(ohms) = %f', Rs(run));
        fprintf('\n Rp(ohms) = %f', Rp(run));
        fprintf('\n RMSE = %f *10^-3\n', RMSE(run)*10^3);
    end
    result(i).data = [Iph, I0, n, Rs, Rp, RMSE];
    result(i).converg = converg;
    result(i).name = selAlgo{i};   
end

%% Tratar cuvas de convergência com tamanhos diferentes
for i = 1: length(result)
    qtdEle = 0;
    for j = 1:length([result(i).converg])
        qtdEle(j) = numel([result(i).converg(j).RMSE]);
    end
    maior = max(qtdEle);
    matrizRMSE = NaN(maior, maior);
    matrizfes = NaN(maior, maior);
    for j = 1:length([result(i).converg])
        matrizRMSE(1:qtdEle(j), j) = [result(i).converg(j).RMSE];
        matrizfes(1:qtdEle(j), j) = [result(i).converg(j).fes];
    end
    result(i).matrizRMSE = matrizRMSE;
    result(i).matrizfes = matrizfes;
end

%% plotar curvas de convergência
figure
for i = 1: length(result)
    converg_curve_mean = mean([result(i).matrizRMSE], 2, 'omitnan');
    fes = mean([result(i).matrizfes], 2, 'omitnan');
    semilogy(fes, converg_curve_mean);
    hold on
end
legend(selAlgo);
xlabel('fes', 'FontSize',18);
ylabel('RMSE(log)', 'FontSize',18);

%% Plotar boxplot do RMSE
figure
for i = 1: length(result)
    values = result(i).data;
    vRMSE(:,i) = values(:,6);
end
    boxplot(vRMSE, selAlgo)

% ylim([0.8, 2.2]*10^-3)

% fprintf(' Parametros MÍNIMOS obtidos:')
% fprintf('\n Iph(A) = %f', min(vIph));
% fprintf('\n I0(uA) = %f', min(vI0)*10^6);
% fprintf('\n n = %f', min(vn));
% fprintf('\n Rs(ohms) = %f', min(vRs));
% fprintf('\n Rp(ohms) = %f', min(vRp));
% fprintf('\n RMSE = %f *10^-3\n', min(vRMSE)*10^3);
% vetorDEminimo = [min(vIph);min(vI0)*10^6;min(vn);min(vRs);min(vRp); min(vRMSE)];
% 
% fprintf(' Parametros MÁXIMOS obtidos:')
% fprintf('\n Iph(A) = %f', max(vIph));
% fprintf('\n I0(uA) = %f', max(vI0)*10^6);
% fprintf('\n n = %f', max(vn));
% fprintf('\n Rs(ohms) = %f', max(vRs));
% fprintf('\n Rp(ohms) = %f', max(vRp));
% fprintf('\n RMSE = %f *10^-3\n', max(vRMSE)*10^3);
% vetorDEmaximo = [max(vIph);max(vI0)*10^6;max(vn);max(vRs);max(vRp); max(vRMSE)];
% 
% fprintf(' Parametros MÉDIOS obtidos:')
% fprintf('\n Iph(A) = %f', mean(vIph));
% fprintf('\n I0(uA) = %f', mean(vI0)*10^6);
% fprintf('\n n = %f', mean(vn));
% fprintf('\n Rs(ohms) = %f', mean(vRs));
% fprintf('\n Rp(ohms) = %f', mean(vRp));
% fprintf('\n RMSE = %f *10^-3\n', mean(vRMSE)*10^3);
% vetorDEmedia = [mean(vIph);mean(vI0)*10^6;mean(vn);mean(vRs);mean(vRp);mean(vRMSE)];
% 
% fprintf('\n Desvio Padrao:')
% fprintf('\n std_Iph(A) = %f', std(vIph));
% fprintf('\n std_I0(uA) = %f', std(vI0)*10^6);
% fprintf('\n std_n = %f', std(vn));
% fprintf('\n std_Rs(ohms) = %f', std(vRs));
% fprintf('\n std_Rp(ohms) = %f', std(vRp));
% fprintf('\n std_RMSE = %f *10^-3\n', std(vRMSE)*10^3);
% vetorDEstd = [std(vIph);std(vI0)*10^6;std(vn);std(vRs);std(vRp); std(vRMSE)];
% vetorDEtudo = [vetorDEminimo, vetorDEmaximo, vetorDEmedia, vetorDEstd];
% 
% fprintf('\n Tempo médio de execução = %f \n', mean(elapseTime));

