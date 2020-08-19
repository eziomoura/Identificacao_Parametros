function [xBest, fBest, fBestCurve, fesCurve] = TLABC(fobj, LB, UB, POP_SIZE, MAX_FES, SHOW_CONVERG)
% Descrição
%   TLABC minimiza a fobj usando a metaheurística "Teaching–learning–based artificial bee colony",
%  conforme descrito em [1]. Autor não indica uma estratégia para tratar os
%  elementos que saiam do espaço de busca.
%
% Entradas:
%   fobj - Função objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada parâmetro
%   UB - Vetor linha com os limites superior de cada parâmetro
%   POP_SIZE - Inteiro com o tamanho da população
%   MAX_FES - Inteiro com o quantidade máxima de avalições da função objetivo
%   SHOW_CONVERG - Valor boleador que se for VERDADEIRO, ativará as saídas com os vetores
%           referentes a curva de convergêngia (converg_RMSE e converg_fes)
%
% Saídas:
%   xBest - Vetor com os parâmetros que minimizam fobj
%   fBest - Valor da fobj avaliada em xBest
%   fBestCurve - Vetor com o fBest ao final de cada iteração
%   fesCurve - Vetor com o número de avalições  da função objetivo ao
%       final de cada iteração
%
% Fontes:
%   [1]CHEN, X.et al. Teaching–learning–based artificial bee colony for solar photovoltaic parameter estimation. Applied Energy, Elsevier, v. 212, n. December 2017, p. 1578–1588, 2018. ISSN 03062619. Disponível em: <https://doi.org/10.1016/j.apenergy.2017.12.115>

%% parâmetros do algoritmo
DIM = length(LB); % qtd de variaveis de design

% numero de tentativas de melhoramento de uma solução candidata
LIMIT = 200;
F = rand; % fator de escala

% inicialização da população
x = LB + (UB - LB).*rand(POP_SIZE, DIM); % População
fobjValue = fobj(x);           % Valor da função objetivo para cada solução candidata
fitValue = fitness(fobjValue); % Fitness de cada x. Quanto maior melhor
trial = zeros(1,POP_SIZE);    % Quantidade de tentativas de melhoramento de uma solução
fes = POP_SIZE;               % Quantidade de avalições da função objetivo
%% pre alocacao da curva de convergência
if SHOW_CONVERG
    %    MAX_ITER = floor((MAX_FES - POP_SIZE)/(2*POP_SIZE));
    %     fBestCurve = zeros(MAX_ITER +1,1);
    %     fesCurve = zeros(MAX_ITER +1,1);
    fBestCurve(1) = min(fobjValue);
    fesCurve(1) = fes;
end
iter = 1; % contador de iterações
while(fes + 2*POP_SIZE + 1 <= MAX_FES)
    %% Fase Teaching-based employed bee
    for i = 1:POP_SIZE
        xMean = mean(x);
        
        % determina o professor
        [~,idBest] = min(fobjValue);
        xBest = x(idBest,:);
        
        % Teaching factor
        TF = randi([1 2]);
        
        % selecionar 3 abelhas diferentes de i
        id = randperm(POP_SIZE,4);
        id(i==id) = [];
        
        % produzir nova solução
        for d = 1:DIM
            if rand < 0.5
                xNew(d) = x(i,d) + rand*(xBest(d) - TF*xMean(d));
            else
                xNew(d) = x(id(1), d) + F*(x(id(2), d) - x(id(3), d));
            end
        end
%         
%         if rand < 0.5
%             xNew = x(i,:) + rand*(xBest - TF*xMean);
%         else
%             xNew = x(id(1), :) + F*(x(id(2),:) - x(id(3),:));
%         end
        
        
        % checar limitantes
        for k = 1:DIM
            if xNew(k) < LB(k)
                xNew(k) = LB(k);
            elseif xNew(k) > UB(k)
                xNew(k) = UB(k);
            end
        end
        
        % avaliar a nova solução
        fobjValueNew = fobj(xNew);
        fitValueNew = fitness(fobjValueNew);
        fes = fes +1;
        
        % quanto maior o fit melhor a solução
        if (fitValueNew > fitValue(i))
            x(i,:) = xNew;
            fitValue(i) = fitValueNew;
            fobjValue(i) = fobjValueNew;
            trial(i) = 0;
        else
            trial(i) = trial(i) + 1;
        end
    end
    
    %% Fase Learning-based on looker bee
    prob = fitValue/sum(fitValue); % probabilidade de uma fonte de comida ser selecionada
    for i = 1:POP_SIZE
        if(prob(i) > rand)
            % Seleciona outra fonte de alimentos
            k = randi(POP_SIZE);
            while(k == i)
                k = randi(POP_SIZE);
            end
            
            % Gera nova solução
            if fobjValue(i) <= fobjValue(k)
                xNew = x(i,:) + rand*(x(i,:) - x(k,:));
            else
                xNew = x(i,:) + rand*(x(k,:) - x(i,:));
            end
            
            % Checa limitantes
            
            for k = 1:DIM
                if xNew(k) < LB(k)
                    xNew(k) = LB(k);
                elseif xNew(k) > UB(k)
                    xNew(k) = UB(k);
                end
            end
            % Avalia nova solução
            fobjValueNew = fobj(xNew);
            fitValueNew = fitness(fobjValueNew);
            fes = fes +1;
            
            % Quanto maior o fit melhor a solução
            if (fitValueNew > fitValue(i))
                x(i,:) = xNew;
                fitValue(i) = fitValueNew;
                fobjValue(i) = fobjValueNew;
                trial(i) = 0;
            else
                trial(i) = trial(i) + 1;
            end
        end
    end
    
    %% SCOUT BEE
    % Há no máximo uma scout por fase
    [~, ind] =  max(trial);
    if (trial(ind) > LIMIT)
        trial(ind) = 0;
        xRand = LB + (UB-LB).*rand(1,DIM);
        xGO = rand*(max(x) + min(x)) - xRand;
        fRand = fobj(xRand);
        fGO = fobj(xGO);
        if fGO < fRand
            x(ind,:) = xGO;
            fobjValue(ind) = fGO;
            fitValue(ind) = fitness(fGO);
        else
            x(ind,:) = xRand;
            fobjValue(ind) = fRand;
            fitValue(ind) = fitness(fRand);
        end
        fes = fes +1;
    end
    
    % Atualiza curva de convergencia
    iter = iter +1;
    if SHOW_CONVERG
        fBestCurve(iter,1) = min(fobjValue);
        fesCurve(iter,1) = fes;
    end
end
[fBest, id] = min(fobjValue);
xBest = x(id,:);
end

function fitval = fitness(fobjValue)
fitval = zeros(size(fobjValue));
id = (fobjValue >= 0);
fitval(id) = 1./(fobjValue(id)+1);
fitval(~id)  = 1 + abs(fobjValue(~id));
end