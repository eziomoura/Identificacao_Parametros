function [xBest, fBest, fBestCurve, fesCurve] = ITLBO(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
% Descrição
%     ITLBO minimiza a fobj usando a metaheurística "Improved Teaching-learning-based",
% conforme descrita em [1]. Como em [1] não é explicitado a forma de 
% tratamento das restrições, aqui será atribuido um valor aleatório ao
% parâmetro que ultrapasse seu limite superior ou inferior.
%
% Entradas:
%   fobj - Função objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada parâmetro
%   UB - Vetor linha com os limites superior de cada parâmetro
%   POP_SIZE - Inteiro com o tamanho da população
%   MAX_FES - Inteiro com o quantidade máxima de avalições da função objetivo
%   SHOW_CONVERG - Valor boleano que se for VERDADEIRO, ativará as saídas com os vetores 
%       referentes a curva de convergêngia (converg_RMSE e converg_fes)
%        
% Saídas:
%   xBest - Vetor com os parâmetros que minimizam fobj
%   fBest - Valor da fobj avaliada em xBest
%   fBestCurve - Vetor com o fBest ao final de cada iteração
%   fesCurve - Vetor com o número de avalições  da função objetivo ao
%       final de cada iteração
%
% Fontes:
%   [1] LI, S.; GONG, W.; YAN, X.; HU, C.; BAI, D.; WANG, L.; GAO, L. Parameter extraction of photovoltaic models using an improved teaching-learning-based optimization. Energy Conversion and Management, v. 186, n. December 2018, p. 293–305, 2019. 
%%
% Parametros do algoritmo
POP_SIZE = PARAM.pop;

% Generation of initial population
DIM = length(LB);     % qtd de variaveis de design
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
fit = fobj(x);

fes = POP_SIZE; % contador de avaliacoes da funcao objetivo

% pre-alocacao de memoria
MAX_ITER = floor((MAX_FES - POP_SIZE)/(2*POP_SIZE));

% pre alocacao da curva de convergência
if SHOW_CONVERG
    fBestCurve = zeros(MAX_ITER + 1, 1);
    fesCurve =  zeros(MAX_ITER + 1, 1);
    fBestCurve(1) = min(fit);
    fesCurve(1) = fes;
end

iter = 1; % contador de iteracoes
while(fes+2*POP_SIZE <= MAX_FES)
    for i = 1:POP_SIZE
        %% Teacher Phase
        xMean = mean(x);
        fitMean = fobj(xMean);
        
        % Determination of teacher
        [~,id] = min(fit);
        xBest = x(id,:);
        
        % Teaching factor
        TF = randi([1 2]);
        
        % Generation of a new solution
        if fit(i) < fitMean
            idRand = randperm(POP_SIZE, 3);
            idRand(idRand == i) = [];
            xNew = x(i,:) + rand*(xBest - x(i,:))...
                          + rand*(x(idRand(1),:) - x(idRand(2),:));
        else
            xNew = x(i,:) + rand*(xBest - TF*xMean);
        end
        
        % Checar limites
        xNew = boudaryCorrection(xNew, LB, UB, DIM, 1);        
        
        % avaliar a nova solucao
        fitNew = fobj(xNew);
        
        % Greedy selection
        if (fitNew < fit(i))
            x(i,:) = xNew;
            fit(i) = fitNew;
        end
%%      Learner Phase
        % selecionar parceiros aleatoriamente
        idRand = randperm(POP_SIZE,5);
        idRand(idRand==i) = [];
        
        % Generation of a new solution
        if (fit(i) < fitMean)
            xNew = x(i,:) + rand*(x(idRand(1),:) - x(idRand(2),:));
        else
            xNew = x(i,:) + rand*(x(idRand(1),:) - x(idRand(2),:))...
                          + rand*(x(idRand(3),:) - x(idRand(4),:));
        end
        
        % Checar limites
        xNew = boudaryCorrection(xNew, LB, UB, DIM, 1);
        
        % Avaliar funçaõ objetivo
        fitNew = fobj(xNew);
        
        % Greedy selection
        if(fitNew <  fit(i))
            x(i,:) = xNew;
            fit(i) = fitNew;
        end        
    end
    
    fes = fes + 2*POP_SIZE;
    iter = iter + 1;
    
    if SHOW_CONVERG
        fBestCurve(iter) = min(fit);
        fesCurve(iter) = fes;
    end  
end
% Melhor solução
[fBest, id] = min(fit);
xBest = x(id,:);
end
%%
function xNew = boudaryCorrection(xNew, LB, UB, DIM, POP_SIZE)
u = (xNew < LB) | (xNew > UB);
randomMatrix = LB + (UB - LB).*rand(POP_SIZE, DIM);
xNew(u) = randomMatrix(u);
end
