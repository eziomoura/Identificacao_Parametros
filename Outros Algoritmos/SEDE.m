function [xBest, fBest, fBestCurve, fesCurve] = SEDE(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
% Descrição
%     SEDE minimiza a fobj usando a metaheurística "self-adaptive ensemble-based
% differential evolution" conforme descrita em [1].
% Em [2] é dito que "If the boundary constraint is violated, set the value to be the middle
% of the previous value and the bound"
% 
% Entradas:
%   fobj - Função objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada parâmetro
%   UB - Vetor linha com os limites superior de cada parâmetro
%   PARAM - Estrutura com o seguintes campos:
%      pop - Tamanho da população
%   MAX_FES - Inteiro com o quantidade máxima de avalições da função objetivo
%   SHOW_CONVERG - Valor boleano que se for VERDADEIRO, ativará as saídas com os vetores 
%       referentes a curva de convergência (converg_RMSE e converg_fes)
%        
% Saídas:
%   xBest - Vetor com os parâmetros que minimizam fobj
%   fBest - Valor da fobj avaliada em xBest
%   fBestCurve - Vetor com o fBest ao final de cada iteração
%   fesCurve - Vetor com o número de avalições  da função objetivo ao
%       final de cada iteração
%
% Fontes:
%   [1] LIANG, J.; QIAO, K.; YU, K.; GE, S.; QU, B.; XU, R.; LI, K. Parameters estimation of solar photovoltaic models via a self-adaptive ensemble-based differential evolution. Solar Energy, v. 207, n. C, p. 336–346, 2020.
%   [2] Source code provided by the authors at: https://github.com/cilabzzu/Publications/tree/master/SEDE
%% parâmetros do algoritmo
POP_SIZE = PARAM.pop;
% três conjuntos de parametros
F  = [1.0 1.0 0.8];
CR = [0.1 0.9 0.2];
%% Inicializa a populacao
DIM = length(LB); % qtd de variaveis de design
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
fit = fobj(x);
fes = POP_SIZE;
iter = 1;

if SHOW_CONVERG
    MAX_ITER = floor((MAX_FES-fes)/POP_SIZE) + 1;
    fBestCurve = NaN(MAX_ITER,1);
    fesCurve = NaN(MAX_ITER,1);
    fBestCurve(1) = min(fit);
    fesCurve(1) = fes;
end

while(fes + POP_SIZE <= MAX_FES)
    idParam = randi(3);
    [~, idBest] = min(fit);
    % Mutação e Crossover
    for i = 1:POP_SIZE
        % seleciona 3 individuos
        id = randperm(POP_SIZE,4);
        id(id == i) = [];
        
        if rand > (fes/MAX_FES)
        %% Group1: exploration group

            if  rand < 0.5
                % DE/rand/1
                xNew = x(id(1),:) + F(idParam)*(x(id(2),:) - x(id(3),:));
                
                % Crossover Binomial
                jRand = randi(DIM);
                for j = 1:DIM
                    if rand < CR(idParam) | j == jRand
                        xNew(j) = xNew(j);
                    else
                        xNew(j) = x(i,j);
                    end
                end
            else
                % "current to rand/1", segundo ref[2]
                xNew = x(i, :) + rand*(x(id(1),:) - x(i,:)) + ...
                    F(idParam).*(x(id(2), :) - x(id(3), :));
                % ref2: "Binomial crossover is not used to generate the trial vector under this
                % condition"
            end
        else
        % Group2: exploitation group
            if rand < 0.5
            % "current-to-best"
                xNew = x(i,:) + F(idParam)*(x(idBest,:) - x(i,:)) + F(idParam)*(x(id(1),:) - x(id(2),:));
                
                % Crossover Binomial
                jRand = randi(DIM);
                for j = 1:DIM
                    if rand < CR(idParam) | j == jRand
                        xNew(j) = xNew(j);
                    else
                        xNew(j) = x(i,j);
                    end
                end
            else
            % "current to rand/1", segundo ref[2]
                xNew = x(i, :) + rand*(x(id(1),:) - x(i,:)) + ...
                    F(idParam).*(x(id(2), :) - x(id(3), :));
                
            % ref2: "Binomial crossover is not used to generate the
            % trial vector under this condition"
            end
        end
        
        % checa limitantes
         xNew = boudaryCorrection(xNew, x(i,:), LB, UB);
         
        % avalia a nova posicao
        fitNew = fobj(xNew);
        if fitNew < fit(i)
            x(i,:) = xNew;
            fit(i) = fitNew;
        end
        fes = fes+1;
    end
    %%
    iter = iter+1;
    if SHOW_CONVERG
        fBestCurve(iter,1) = min(fit);
        fesCurve(iter,1) = fes;
    end
end
% Remove NaNs
fBestCurve = fBestCurve(1:iter, 1);
fesCurve = fesCurve(1:iter, 1);

[fBest, id] = min(fit);
xBest = x(id,:);
end
function xNew = boudaryCorrection(xNew, x, LB, UB)
u = (xNew < LB);
MatrizMedia = (LB + x)/2;
xNew(u) = MatrizMedia(u);

u = (xNew > UB);
MatrizMedia = (UB + x)/2;
xNew(u) = MatrizMedia(u);
end