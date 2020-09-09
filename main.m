% OBS2: CIABC foi perda de tempo. Ver notas no arquivo
% OBS: RMSE deve ser retornado como vetor coluna
% OBS3: pode mudar limites para curvas diferentes (photowatt vs stm6-40)
%%
% allOutputs guarda todos os dados necessário à análise
clc; clear; close all;
rng('shuffle');% avoid repeating the same random number arrays when MATLAB restarts
addpath('.\Outros Algoritmos')
addpath('.\Funções Objetivo')
load 'data/allCurves.mat'
%% options
selectedCurves = [1,3];
objetivo.metricas = {'RMSE'};
objetivo.grandezas = {'I',};   % I - current, P- Power, V - Voltage
objetivo.modelos = {'1D','2D'};% 1D - um diodo; 2D - dois diodos      

selAlgo = {'BFS','EJADE','SEDE','PGJAYA', 'ITLBO'};  % Vetor com os algoritmos que deseja avaliar
RUNS = 10;                  % quantidade de execuções distintas
MAX_FES = 100;              % numero maximo de avalicoes da funcao objetivo
paramData;                  % carrega parametros configurados para cada algoritmo

listAlgo = {'BFS','ABC','DE','EJADE','IJAYA','ITLBO','JADE','PGJAYA','PSO','TLBO'}; % (nao atualizada) Lista de todos algoritmos disponíveis

if strcmp(selAlgo{1}, 'all')
    selAlgo = listAlgo;
end

%%
% Lista de modelos de funções objetivo a serem testadas
numObjs = 0;
for i = 1:numel(objetivo.metricas)
    for j = 1:numel(objetivo.grandezas)
        for k = 1:numel(objetivo.modelos)
            numObjs = numObjs + 1;
            obj_code(numObjs).metrica = objetivo.metricas{i};
            obj_code(numObjs).grandeza = objetivo.grandezas{j};
            obj_code(numObjs).modelo = objetivo.modelos{k}; 
        end
    end
end

%% pre alocação
%converg_curve = zeros(RUNS, maxIter);
elapsedTime = zeros(1,RUNS);
Iph = zeros(RUNS,1); I0 = zeros(RUNS,1); n = zeros(RUNS,1);
Rs = zeros(RUNS,1);  Rp = zeros(RUNS,1); f = zeros(RUNS,1);

%% Testes
for i = 1:length(selectedCurves)
    % carrega dados da curva
    Vmed = [IVCurves(selectedCurves(i)).V];   % Vetor de tensoes medidas  [V]
    Imed = [IVCurves(selectedCurves(i)).I];   % Vetor de correntes medidas [A]
    Tc   = [IVCurves(selectedCurves(i)).T];   % Temperatura [ºC]
    Ns   = [IVCurves(selectedCurves(i)).Ns];  % Numero de celulas em serie
    
    % itera sobre as funcões objetivo desejadas
    for numFun = 1:numObjs
        
        % Define a função objetivo
        fun = makeFobj(Vmed, Imed, Ns, Tc, obj_code(numFun));
        fobj = @fun.Objective;
        
        % Define Limites físicos do modelo
        [limite_inf, limite_sup] = getLimit(obj_code(numFun).modelo, Ns);
        
        % itera sobre todos as metaheuristicas
        data = [];
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
                
                if obj_code(numFun).modelo == '1D'
                    Iph(run) = x(1);
                    I0(run) = x(2);
                    n(run) = x(3);
                    Rs(run) = x(4);
                    Rp(run) = x(5);
                elseif obj_code(numFun).modelo == '2D'
                    Iph(run) = x(1);
                    I01(run) = x(2);
                    I02(run) = x(3);
                    n1(run) = x(4);
                    n2(run) = x(5);
                    Rs(run) = x(6);
                    Rp(run) = x(7);
                else 
                    error('verifique o modelo selecionado');
                end
                
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
    allOutputs(i).name = IVCurves(i).name;
    allOutputs(i).model =  model;
end

%% Tratamento dos resultados
% Tratar curvas de convergência com tamanhos diferentes
% adicionando NaN para vetores possuirem mesma dimensão
for i = 1: length(selectedCurves)
    for numFun = 1:numObjs
        result = allOutputs(i).model(numFun).result;
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
        allOutputs(i).model(numFun).result = result;
    end
end

%% Graficos, tabelas e testes estatisticos
iter = 0;
for numFun = 1:numObjs
    tabelaBest = [];
    for iCurve = 1: length(selectedCurves)
        iter = iter + 1;
        result   = allOutputs(iCurve).model(numFun).result;
        modelo   = allOutputs(iCurve).model(numFun).modelo;
        metrica  = allOutputs(iCurve).model(numFun).metrica;
        grandeza = allOutputs(iCurve).model(numFun).grandeza;
        
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
        tabela_melhores.Properties.Description = sprintf('%s %s %s - %s', modelo, metrica, grandeza, IVCurves(selectedCurves(iCurve)).name);
        vetor_de_tabelas_melhores{iter}.tabela = tabela_melhores;
        
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
        tabela_max_min_mean_sd_cpu.Properties.Description = sprintf('%s %s %s - %s', modelo, metrica, grandeza, allOutputs(iCurve).name);
        vetor_de_tabela_max_min_mean_sd_cpu{iter}.tabela = tabela_max_min_mean_sd_cpu;
        
        
        % Teste de wilcoxon (identificar se há diferenças significativas
        % entres os algoritmos de teste
        % H0: BFS não é diferente do outro algoritmo [median(BFS) = median(outro algorimo)]
        % Ha: BFS não é diferente
        % nivel de significancia: 0.05
        alfa = 0.05;
        ref = 'BFS';
        for i = 1:length(result)
            if strcmp({result(i).name}, ref)
                idRef = i;
                break;
            end
        end
        itemp = 0;
        for i = 1:length(result)
            if idRef ~= i
                itemp = itemp + 1;
                p_value(itemp,1) = signrank(result(idRef).f, result(itemp).f);
                if p_value(itemp,1) <= alfa
                    if (median(result(idRef).f - result(itemp).f) > 0)
                        win(itemp,1) = '+';
                    else
                        win(itemp,1) = '-';
                    end
                else
                    win(itemp,1) = '=';
                end
            end
        end
        tabela_wilcoxon = table(p_value, win);
        ref_vs = {result([1:idRef-1, idRef+1:end]).name}.';
        tabela_wilcoxon = addvars(tabela_wilcoxon, ref_vs, 'Before','p_value');
        tabela_wilcoxon.Properties.Description = sprintf('%s %s %s - %s', modelo, metrica, grandeza, allOutputs(iCurve).name);
        vetor_de_tabelas_wilcoxon{iter}.tabela = tabela_wilcoxon;
        

        % Plotar curvas de convergência
        figure
        subplot(2,1,1)
        for iAlgo = 1: length(result)
            converg_curve_mean = mean([result(iAlgo).matrizRMSE], 2, 'omitnan');
            fes = mean([result(iAlgo).matrizfes], 2, 'omitnan');
            semilogy(fes, converg_curve_mean);
            hold on
        end
        legend(selAlgo);
        xlabel('fes', 'FontSize',18);
        yylabel = sprintf('%s(log)', metrica);
        ylabel(yylabel, 'FontSize',18);
        title(sprintf('%s - %s da %s - Curvas de convergência - Modelo: %s', allOutputs(iCurve).name, metrica, grandeza, modelo))
        
        % Plotar boxplot do f
        subplot(2,1,2)
        for iAlgo = 1: length(result)
            vRMSE(:,iAlgo) =  result(iAlgo).f;
        end
        boxplot(vRMSE, selAlgo)
        title(sprintf('%s - %s da %s - boxplot - Modelo: %s', allOutputs(iCurve).name,  metrica, grandeza, modelo))
    end
end
%% Exportar tabelas para excel
% limpa tabelas
delete('.\resultados\melhores.xlsx','.\resultados\minmaxsdtime.xlsx','.\resultados\wilcoxon_test.xlsx');
% tabela melhores
for i = 1:length(vetor_de_tabelas_melhores)
    writetable(vetor_de_tabelas_melhores{i}.tabela,'.\resultados\melhores.xlsx','Sheet',vetor_de_tabelas_melhores{i}.tabela.Properties.Description);
end
% tabela min max sd time
for i = 1:length(vetor_de_tabela_max_min_mean_sd_cpu)
    writetable(vetor_de_tabela_max_min_mean_sd_cpu{i}.tabela,'.\resultados\minmaxsdtime.xlsx','Sheet',vetor_de_tabela_max_min_mean_sd_cpu{i}.tabela.Properties.Description);
end
% tabela wilcoxon
for i = 1:length(vetor_de_tabelas_wilcoxon)
    writetable(vetor_de_tabelas_wilcoxon{i}.tabela,'.\resultados\wilcoxon_test.xlsx','Sheet',vetor_de_tabelas_wilcoxon{i}.tabela.Properties.Description);
end