%%
% Arquivo demonstracao com codigo que testa os diversos algoritmos
% implementados aqui
%
%
% OBS: RMSE deve ser retornado como vetor coluna
%%
% allOutputs guarda todos os dados necessário à análise
clc; clear; close all;
rng('shuffle');% avoid repeating the same random number arrays when MATLAB restarts
addpath('.\Outros Algoritmos')
addpath('.\Funções Objetivo')
load 'data/todasCurvassel.mat'
%% options
selectedCurves = 1:length(IVCurves);
objetivo.metricas = {'RMSE'};
objetivo.grandezas = {'I',};   % I - current, P- Power, V - Voltage
objetivo.modelos = {'1D'};     % 1D - um diodo; 2D - dois diodos      

selAlgo = {'all'};
%selAlgo = {'BFS','MADE', 'SEDE', 'ITLBO'};% Vetor com os algoritmos que deseja avaliar
RUNS = 30;                  % quantidade de execuções distintas
MAX_FES = 80e3;             % numero maximo de avalicoes da funcao objetivo
paramData;                  % carrega parametros configurados para cada algoritmo

listAlgo = {'BFS','SHADE', 'MADE', 'SEDE', 'EJADE', 'TLBO', 'ITLBO', 'TLABC', 'ABC', 'CIABC', 'PSO', 'ELPSO', 'IJAYA', 'PGJAYA'}; % (nao atualizada) Lista de todos algoritmos disponíveis
%listAlgoTOTAL = {'BFS','SHADE', 'MADE', 'SEDE', 'EJADE', 'DE', 'TLBO', 'ITLBO', 'TLABC', 'ABC', 'CIABC', 'PSO', 'ELPSO', 'IJAYA', 'PGJAYA'}; % (nao atualizada) Lista de todos algoritmos disponíveis
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

%% Testes para cada curva
for i = 1:length(selectedCurves)
    % carrega dados da curva
    Vmed = [IVCurves(selectedCurves(i)).V];   % Vetor de tensoes medidas  [V]
    Imed = [IVCurves(selectedCurves(i)).I];   % Vetor de correntes medidas [A]
    Tc   = [IVCurves(selectedCurves(i)).T];   % Temperatura [ºC]
    Ns   = [IVCurves(selectedCurves(i)).Ns];  % Numero de celulas em serie
    nome_dispositivo = [IVCurves(selectedCurves(i)).name];
    fprintf('Testes para a curva %s', nome_dispositivo);
    
    % itera sobre as funcões objetivo desejadas
    for numFun = 1:numObjs
        
        % Define a função objetivo
        fun = makeFobj(Vmed, Imed, Ns, Tc, obj_code(numFun));
        fobj = @fun.Objective;
        
        % Define Limites físicos do modelo
        [limite_inf, limite_sup] = getLimit(obj_code(numFun).modelo, Ns,IVCurves(selectedCurves(i)));
        
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
                    
                    fprintf(' Iph(A) = %f', Iph(run));
                    fprintf('\n I0(uA) = %f', I0(run)*10^6);
                    fprintf('\n n = %f', n(run));
                    fprintf('\n Rs(ohms) = %f', Rs(run));
                    fprintf('\n Rp(ohms) = %f', Rp(run));
                elseif obj_code(numFun).modelo == '2D'
                    Iph(run) = x(1);
                    I01(run) = x(2);
                    I02(run) = x(3);
                    n1(run) = x(4);
                    n2(run) = x(5);
                    Rs(run) = x(6);
                    Rp(run) = x(7);
                    
                     
                    fprintf(' Iph(A) = %f', Iph(run));
                    fprintf('\n I01(uA) = %f', I01(run)*10^6);
                    fprintf('\n I02(uA) = %f', I02(run)*10^6);
                    fprintf('\n n1 = %f', n1(run));
                    fprintf('\n n2 = %f', n2(run));
                    fprintf('\n Rs(ohms) = %f', Rs(run));
                    fprintf('\n Rp(ohms) = %f', Rp(run));
                else 
                    error('verifique o modelo selecionado');
                end
                fprintf('\n RMSE = %f *10^-3\n', f(run)*10^3);
                fprintf('tempo decorrido na estimacao %d', elapsedTime(run));
                
            end
            % result guarda o resultado de todas metaheuristica para uma
            % curva especifica com modelo especificado
            result_metaheuristc(iAlgo).x = data;
            result_metaheuristc(iAlgo).f = f;
            result_metaheuristc(iAlgo).duration = elapsedTime;
            result_metaheuristc(iAlgo).converg = converg;
            result_metaheuristc(iAlgo).name = selAlgo{iAlgo};
        end % fim loop por todas metaheuristicas
        % model guarda o resultado de todas metaheuristica para uma
        % curva especifica com modelo especificado
        model(numFun).result = result_metaheuristc;
        model(numFun).metrica = fun.metrica;
        model(numFun).grandeza = fun.grandeza;
        model(numFun).modelo = fun.modelo;
    end % fim loop por todos modelos
    allOutputs(i).name = IVCurves(selectedCurves(i)).name;
    allOutputs(i).model =  model;
end% fim loop por todas as curvas

%% Tratamento dos resultados
% Tratar curvas de convergência com tamanhos diferentes
% adicionando NaN para vetores possuirem mesma dimensão
for i = 1: length(selectedCurves)
    for numFun = 1:numObjs
        result_metaheuristc = allOutputs(i).model(numFun).result;
        for iAlgo = 1: length(result_metaheuristc)
            qtdEle = 0;
            for j = 1:length([result_metaheuristc(iAlgo).converg])
                qtdEle(j) = numel([result_metaheuristc(iAlgo).converg(j).f]);
            end
            maior = max(qtdEle);
            matrizRMSE = NaN(maior, maior);
            matrizfes = NaN(maior, maior);
            for j = 1:length([result_metaheuristc(iAlgo).converg])
                matrizRMSE(1:qtdEle(j), j) = [result_metaheuristc(iAlgo).converg(j).f];
                matrizfes(1:qtdEle(j), j)  = [result_metaheuristc(iAlgo).converg(j).fes];
            end
            result_metaheuristc(iAlgo).matrizRMSE = matrizRMSE;
            result_metaheuristc(iAlgo).matrizfes = matrizfes;
        end
        allOutputs(i).model(numFun).result = result_metaheuristc;
    end
end

%% Realiza testes estatisticos e gera tabelas com resultados
iter = 0;
for numFun = 1:numObjs
    tabelaBest = [];
    for iCurve = 1: length(selectedCurves)
        iter = iter + 1;
        result_metaheuristc   = allOutputs(iCurve).model(numFun).result;
        modelo   = allOutputs(iCurve).model(numFun).modelo;
        metrica  = allOutputs(iCurve).model(numFun).metrica;
        grandeza = allOutputs(iCurve).model(numFun).grandeza;
        
        % Tabela com melhor f(x_best) de cada metaheuristica
        for iAlgo = 1: length(result_metaheuristc)
            [fval, idBest] = min(result_metaheuristc(iAlgo).f);
            tabelaBest(iAlgo,:) = [result_metaheuristc(iAlgo).x(idBest,:)];
            fbest(iAlgo,1) = fval;
        end
        if strcmp(modelo, '1D')
            labels = {'Iph', 'I0', 'n', 'Rs', 'Rp'};
        else
            labels = {'Iph', 'I01', 'I02', 'n1','n2', 'Rs', 'Rp'};
        end
        tabela_melhores = array2table(tabelaBest, 'VariableNames', labels);
        Algorithm = {result_metaheuristc(:).name}.';
        tabela_melhores = addvars(tabela_melhores, Algorithm, 'Before','Iph');
        tabela_melhores = addvars(tabela_melhores, fbest, 'After','Rp', 'NewVariableNames', metrica);
        tabela_melhores.Properties.Description = sprintf('%s %s %s - %s', modelo, metrica, grandeza, IVCurves(selectedCurves(iCurve)).name);
        vetor_de_tabelas_melhores{iter}.tabela = tabela_melhores;
        
        % Tabela max, min, mean, sd and CPUtime
        for iAlgo = 1: length(result_metaheuristc)
            fmin(iAlgo,1)  = min(result_metaheuristc(iAlgo).f);
            fmax(iAlgo,1)  = max(result_metaheuristc(iAlgo).f);
            fmean(iAlgo,1) = mean(result_metaheuristc(iAlgo).f);
            fstd(iAlgo,1)  = std(result_metaheuristc(iAlgo).f);
            ftime(iAlgo,1) = mean(result_metaheuristc(iAlgo).duration);
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
        p_value = []; 
        for i = 1:length(result_metaheuristc)
            if strcmp({result_metaheuristc(i).name}, ref)
                idRef = i;
                break;
            end
        end
        itemp = 1;
        for i = 1:length(result_metaheuristc)
            if  i ~= idRef
                p_value(itemp,1) = ranksum(result_metaheuristc(idRef).f, result_metaheuristc(i).f);
                if p_value(itemp,1) <= alfa
                    if (median(result_metaheuristc(idRef).f - result_metaheuristc(i).f) > 0)
                        win(itemp,1) = '-';
                    else
                        win(itemp,1) = '+';
                    end
                else
                    win(itemp,1) = '=';
                end
                itemp = itemp+1;
            end
        end
        tabela_wilcoxon = table(p_value, win);
        ref_vs = {result_metaheuristc([1:idRef-1, idRef+1:end]).name}.';
        tabela_wilcoxon = addvars(tabela_wilcoxon, ref_vs, 'Before','p_value');
        tabela_wilcoxon.Properties.Description = sprintf('%s %s %s - %s', modelo, metrica, grandeza, allOutputs(iCurve).name);
        vetor_de_tabelas_wilcoxon{iter}.tabela = tabela_wilcoxon;
        
       end
end
%% Exportar tabelas para excel
% limpa tabelas
%delete('.\resultados\melhores.xlsx','.\resultados\minmaxsdtime.xlsx','.\resultados\wilcoxon_test.xlsx');

destdirectory = [pwd, '\resultados\', strrep(datestr(datetime), ':', '_')];
mkdir(destdirectory);   %create the directory

% tabela melhores
nomearquivo = 'melhores.xlsx';
fulldestination = fullfile(destdirectory, nomearquivo);
for i = 1:length(vetor_de_tabelas_melhores)
    writetable(vetor_de_tabelas_melhores{i}.tabela,fulldestination,'Sheet',vetor_de_tabelas_melhores{i}.tabela.Properties.Description);
end

% tabela min max sd time
nomearquivo = 'minmaxsdtime.xlsx';
fulldestination = fullfile(destdirectory, nomearquivo);
for i = 1:length(vetor_de_tabela_max_min_mean_sd_cpu)
    writetable(vetor_de_tabela_max_min_mean_sd_cpu{i}.tabela,fulldestination,'Sheet',vetor_de_tabela_max_min_mean_sd_cpu{i}.tabela.Properties.Description);
end

% tabela wilcoxon
nomearquivo = 'wilcoxon_test.xlsx';
fulldestination = fullfile(destdirectory, nomearquivo);
for i = 1:length(vetor_de_tabelas_wilcoxon)
    writetable(vetor_de_tabelas_wilcoxon{i}.tabela,fulldestination,'Sheet',vetor_de_tabelas_wilcoxon{i}.tabela.Properties.Description);
end
% todos os dados
nomearquivo = 'OUTPUT.mat';
fulldestination = fullfile(destdirectory, nomearquivo);
save(fulldestination, '-v7.3')
%plotagens