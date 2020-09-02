% OBS2: CIABC foi perda de tempo. Ver notas no arquivo
% OBS: RMSE deve ser retornado como vetor coluna

clc; clear; close all;
rng('shuffle');% avoid repeating the same random number arrays when MATLAB restarts
addpath('.\Outros Algoritmos')
addpath('.\Funções Objetivo')
load 'data/allCurves.mat'
%%
selectedCurves = [1,2];
CODE_FUN_OBJ = [1,2]; % ver arquivo makeFunObj, para lista de codigos
selAlgo = { 'BFS','EJADE'}; % Vetor com os algoritmos que deseja avaliar
RUNS = 10; % quantidade de execuções distintas
MAX_FES = 100; % numero maximo de avalicoes da funcao objetivo
paramData; % carrega parametros configurados para cada algoritmo
% graphic = false; % deseja plotar curvas IV?
% plotCurvesEveryRun = false;

listAlgo = {'BFS','ABC','DE','EJADE','IJAYA','ITLBO','JADE','PGJAYA','PSO','TLBO'}; % (nao atualizada) Lista de todos algoritmos disponíveis
if strcmp(selAlgo{1}, 'all')
    selAlgo = listAlgo;
end


%% pre alocação
%converg_curve = zeros(RUNS, maxIter);
elapsedTime = zeros(1,RUNS);
Iph = zeros(RUNS,1); I0 = zeros(RUNS,1); n = zeros(RUNS,1);
Rs = zeros(RUNS,1);  Rp = zeros(RUNS,1); RMSE = zeros(RUNS,1);

%% Dados de entrada
for i = 1: length(selectedCurves)
    % carrega dados da curva
    Vmed = [IVCurve(i).V];   % Vetor de tensoes medidas  [V]
    Imed = [IVCurve(i).I];   % Vetor de correntes medidas [A]
    Tc   = [IVCurve(i).T];   % Temperatura [ºC]
    Ns   = [IVCurve(i).Ns];  % Numero de celulas em serie
    
    % itera sobre as funcões objetivo desejadas
    for numFun = 1:length(CODE_FUN_OBJ)
        
        % Define a função objetivo
        fun = makeFobj(Vmed, Imed, Ns, Tc, CODE_FUN_OBJ(numFun));
        fobj = @fun.Objective;
        
        % Define Limites físicos do modelo
        [limite_inf, limite_sup] = getLimit(funType, Ns);
        
        % itera sobre todos as metaheuristicas
        for iAlgo = 1: length(selAlgo)
            fprintf('\nTestando %s \n', selAlgo{iAlgo});
            
            for run = 1:RUNS
                [metaheuristic, prmt] = getAlgo(selAlgo{iAlgo}, Param);
                fprintf('\nExecução %d \n', run)
                tic
                [x, RMSE(run), converg(run).RMSE, converg(run).fes] =  metaheuristic(fobj, limite_inf, limite_sup,...
                    prmt, MAX_FES, true);
                
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
            result(iAlgo).data = [Iph, I0, n, Rs, Rp, RMSE];
            result(iAlgo).converg = converg;
            result(iAlgo).name = selAlgo{iAlgo};
        end
        model(numFun).result = result;
        model(numFun).Fobj = CODE_FUN_OBJ(numFun);
    end
end



%% Tratamento dos resultados
% Tratar curvas de convergência com tamanhos diferentes
for iAlgo = 1: length(result)
    qtdEle = 0;
    for j = 1:length([result(iAlgo).converg])
        qtdEle(j) = numel([result(iAlgo).converg(j).RMSE]);
    end
    maior = max(qtdEle);
    matrizRMSE = NaN(maior, maior);
    matrizfes = NaN(maior, maior);
    for j = 1:length([result(iAlgo).converg])
        matrizRMSE(1:qtdEle(j), j) = [result(iAlgo).converg(j).RMSE];
        matrizfes(1:qtdEle(j), j) = [result(iAlgo).converg(j).fes];
    end
    result(iAlgo).matrizRMSE = matrizRMSE;
    result(iAlgo).matrizfes = matrizfes;
end

% plotar curvas de convergência
figure
for iAlgo = 1: length(result)
    converg_curve_mean = mean([result(iAlgo).matrizRMSE], 2, 'omitnan');
    fes = mean([result(iAlgo).matrizfes], 2, 'omitnan');
    semilogy(fes, converg_curve_mean);
    hold on
end
legend(selAlgo);
xlabel('fes', 'FontSize',18);
ylabel('RMSE(log)', 'FontSize',18);

% Plotar boxplot do RMSE
figure
for iAlgo = 1: length(result)
    values = result(iAlgo).data;
    vRMSE(:,iAlgo) = values(:,6);
end
boxplot(vRMSE, selAlgo)






%% etc
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

