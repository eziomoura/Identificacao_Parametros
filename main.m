% OBS2: CIABC foi perda de tempo. Ver notas no arquivo
% OBS: RMSE deve ser retornado como vetor coluna

clc; clear; close all;
rng('shuffle');% avoid repeating the same random number arrays when MATLAB restarts
addpath('.\Outros Algoritmos')
addpath('.\Funções Objetivo')
load 'data/allCurves.mat'
%% options
selectedCurves = [1,2];
CODE_FUN_OBJ = [1,2];       % ver arquivo makeFunObj, para lista de codigos
selAlgo = {'BFS','EJADE'};  % Vetor com os algoritmos que deseja avaliar
RUNS = 10;                  % quantidade de execuções distintas
MAX_FES = 100;            % numero maximo de avalicoes da funcao objetivo
paramData;                  % carrega parametros configurados para cada algoritmo

listAlgo = {'BFS','ABC','DE','EJADE','IJAYA','ITLBO','JADE','PGJAYA','PSO','TLBO'}; % (nao atualizada) Lista de todos algoritmos disponíveis
if strcmp(selAlgo{1}, 'all')
    selAlgo = listAlgo;
end


%% pre alocação
%converg_curve = zeros(RUNS, maxIter);
elapsedTime = zeros(1,RUNS);
Iph = zeros(RUNS,1); I0 = zeros(RUNS,1); n = zeros(RUNS,1);
Rs = zeros(RUNS,1);  Rp = zeros(RUNS,1); f = zeros(RUNS,1);

%% Testes
for i = 1: length(selectedCurves)
    % carrega dados da curva
    Vmed = [IVCurves(i).V];   % Vetor de tensoes medidas  [V]
    Imed = [IVCurves(i).I];   % Vetor de correntes medidas [A]
    Tc   = [IVCurves(i).T];   % Temperatura [ºC]
    Ns   = [IVCurves(i).Ns];  % Numero de celulas em serie
    
    % itera sobre as funcões objetivo desejadas
    for numFun = 1:length(CODE_FUN_OBJ)
        
        % Define a função objetivo
        fun = makeFobj(Vmed, Imed, Ns, Tc, CODE_FUN_OBJ(numFun));
        fobj = @fun.Objective;
        
        % Define Limites físicos do modelo
        [limite_inf, limite_sup] = getLimit(CODE_FUN_OBJ(numFun), Ns);
        
        % itera sobre todos as metaheuristicas
        for iAlgo = 1: length(selAlgo)
            fprintf('\nTestando %s \n', selAlgo{iAlgo});
            
            for run = 1:RUNS
                [metaheuristic, prmt] = getAlgo(selAlgo{iAlgo}, Param);
                fprintf('\nExecução %d \n', run)
                tic
                [x, f(run), converg(run).f, converg(run).fes] =  metaheuristic(fobj, limite_inf, limite_sup,...
                                                                                       prmt, MAX_FES, true);
                
                elapsedTime(run) = toc;
                data(run,:) = x;
                
                Iph(run) = x(1);
                I0(run) = x(2);
                n(run) = x(3);
                Rs(run) = x(4);
                Rp(run) = x(5);
                fprintf(' Iph(A) = %f', Iph(run));
                fprintf('\n I0(uA) = %f', I0(run)*10^6);
                fprintf('\n n = %f', n(run));
                fprintf('\n Rs(ohms) = %f', Rs(run));
                fprintf('\n Rp(ohms) = %f', Rp(run));
                fprintf('\n RMSE = %f *10^-3\n', f(run)*10^3);
            end
            result(iAlgo).x = data;
            result(iAlgo).f = f;
            result(iAlgo).duration = elapsedTime;
            result(iAlgo).converg = converg;
            result(iAlgo).name = selAlgo{iAlgo};
        end
        model(numFun).result = result;
        model(numFun).metrica = fun.metrica;
        model(numFun).grandeza = fun.grandeza;
        model(numFun).modelo = fun.modelo;
    end
    Benchmark(i).name = IVCurves(i).name;
    Benchmark(i).model =  model;
end



%% Tratamento dos resultados
% Tratar curvas de convergência com tamanhos diferentes
% adicionando NaN para vetores possuirem mesma dimensão
for i = 1: length(selectedCurves)
    for numFun = 1:length(CODE_FUN_OBJ)
        result = Benchmark(i).model(numFun).result;
        for iAlgo = 1: length(result)
            qtdEle = 0;
            for j = 1:length([result(iAlgo).converg])
                qtdEle(j) = numel([result(iAlgo).converg(j).f]);
            end
            maior = max(qtdEle);
            matrizRMSE = NaN(maior, maior);
            matrizfes = NaN(maior, maior);
            for j = 1:length([result(iAlgo).converg])
                matrizRMSE(1:qtdEle(j), j) = [result(iAlgo).converg(j).f];
                matrizfes(1:qtdEle(j), j)  = [result(iAlgo).converg(j).fes];
            end
            result(iAlgo).matrizRMSE = matrizRMSE;
            result(iAlgo).matrizfes = matrizfes;
        end
        Benchmark(i).model(numFun).result = result;
    end
end

%% Graficos, tabelas e testes estatisticos
iter = 1;
for i = 1: length(selectedCurves)
    for numFun = 1:length(CODE_FUN_OBJ)
        result = Benchmark(i).model(numFun).result;
        modelo = Benchmark(i).model(numFun).modelo;
        metrica = Benchmark(i).model(numFun).metrica;
        grandeza = Benchmark(i).model(numFun).grandeza;
        
        % Tabela com melhor f de cada metaheuristica
        for iAlgo = 1: length(result)
            [fval, idBest] = min(result(iAlgo).f);
            tabelaBest(iAlgo,:) = [result(iAlgo).x(idBest,:)];
            fbest(iAlgo,1) = fval;
        end
        if strcmp(modelo, '1D')
            labels = {'Iph', 'I0', 'n', 'Rs', 'Rp'};
        else
            labels = {'Iph', 'I01', 'I02', 'n1','n2', 'Rs', 'Rp'};
        end
        tabela_melhores = array2table(tabelaBest, 'VariableNames', labels);
        Algorithm = {result(:).name}.';
        tabela_melhores = addvars(tabela_melhores, Algorithm, 'Before','Iph');
        tabela_melhores = addvars(tabela_melhores, fbest, 'After','Rp', 'NewVariableNames', metrica);
        tabela_melhores.Properties.Description = sprintf('comparativo - %s %s %s - %s', modelo, metrica, grandeza, IVCurves(selectedCurves(i)).name);
        vetor_de_tabelas_melhores{iter} = tabela_melhores;
        
        % Tabela max, min, mean, sd and CPUtime
        for iAlgo = 1: length(result)
            fmin(iAlgo,1)  = min(result(iAlgo).f);
            fmax(iAlgo,1)  = max(result(iAlgo).f);
            fmean(iAlgo,1) = mean(result(iAlgo).f);
            fstd(iAlgo,1)  = std(result(iAlgo).f);
            ftime(iAlgo,1) = mean(result(iAlgo).duration);
        end
        labels = {'Algorithm', 'Min', 'Max', 'Mean', 'Std', 'CPU time'};
        tabela_max_min_mean_sd_cpu = table(Algorithm, fmin, fmax, fmean, fstd, ftime, 'VariableNames', labels);
        vetor_de_tabela_max_min_mean_sd_cpu{iter} = tabela_max_min_mean_sd_cpu;
        
        
        % Teste de wilcoxon (identificar se há diferenças significativas
        % entres os algoritmos de teste
        % H0: BFS não é diferente do outro algoritmo [median(BFS) = median(outro algorimo)]
        % Ha: BFS não é diferente
        % nivel de significancia: 0.05
        alfa = 0.05;
        idBFS = sum(strcmp({result(:).name}, 'BFS'));
        for i = 1:length(result)
            if i ~= idBFS
                pval = signrank(result(idBFS).f, result(i).f);
                if pval < alfa
                    if median(result(idBFS).f) < median(result(i).f)
                        wilcoxon(i) = '+';
                    else
                        wilcoxon(i) = '-';
                    end
                else
                    wilcoxon(i) = '=';
                end
            end
        end
        

        % Plotar curvas de convergência
        figure
        for iAlgo = 1: length(result)
            converg_curve_mean = mean([result(iAlgo).matrizRMSE], 2, 'omitnan');
            fes = mean([result(iAlgo).matrizfes], 2, 'omitnan');
            semilogy(fes, converg_curve_mean);
            hold on
        end
        legend(selAlgo);
        xlabel('fes', 'FontSize',18);
        ylabel('f(log)', 'FontSize',18);
        title(sprintf('%s - %s da %s - Curvas de convergência - Modelo: %s', IVCurves(selectedCurves(i)).name, metrica, grandeza, modelo))
        
        % Plotar boxplot do f
        figure
        for iAlgo = 1: length(result)
            vRMSE(:,iAlgo) =  result(iAlgo).f;
        end
        boxplot(vRMSE, selAlgo)
        title(sprintf('%s - %s da %s - boxplot - Modelo: %s', IVCurves(selectedCurves(i)).name,  metrica, grandeza, modelo))
    end
end




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

