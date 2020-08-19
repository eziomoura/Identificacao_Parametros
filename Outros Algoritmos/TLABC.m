function [xBest, fBest, fBestCurve, fesCurve] = TLABC(fobj, LB, UB, POP_SIZE, MAX_FES, SHOW_CONVERG)
% Descri��o
%   TLABC minimiza a fobj usando a metaheur�stica "Teaching�learning�based artificial bee colony",
%  conforme descrito em [1]. Autor n�o indica uma estrat�gia para tratar os
%  elementos que saiam do espa�o de busca.
%
% Entradas:
%   fobj - Fun��o objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada par�metro
%   UB - Vetor linha com os limites superior de cada par�metro
%   POP_SIZE - Inteiro com o tamanho da popula��o
%   MAX_FES - Inteiro com o quantidade m�xima de avali��es da fun��o objetivo
%   SHOW_CONVERG - Valor boleador que se for VERDADEIRO, ativar� as sa�das com os vetores
%           referentes a curva de converg�ngia (converg_RMSE e converg_fes)
%
% Sa�das:
%   xBest - Vetor com os par�metros que minimizam fobj
%   fBest - Valor da fobj avaliada em xBest
%   fBestCurve - Vetor com o fBest ao final de cada itera��o
%   fesCurve - Vetor com o n�mero de avali��es  da fun��o objetivo ao
%       final de cada itera��o
%
% Fontes:
%   [1]CHEN, X.et al. Teaching�learning�based artificial bee colony for solar photovoltaic parameter estimation. Applied Energy, Elsevier, v. 212, n. December 2017, p. 1578�1588, 2018. ISSN 03062619. Dispon�vel em: <https://doi.org/10.1016/j.apenergy.2017.12.115>

%% par�metros do algoritmo
DIM = length(LB); % qtd de variaveis de design

% numero de tentativas de melhoramento de uma solu��o candidata
LIMIT = 200;
F = rand; % fator de escala

% inicializa��o da popula��o
x = LB + (UB - LB).*rand(POP_SIZE, DIM); % Popula��o
fobjValue = fobj(x);           % Valor da fun��o objetivo para cada solu��o candidata
fitValue = fitness(fobjValue); % Fitness de cada x. Quanto maior melhor
trial = zeros(1,POP_SIZE);    % Quantidade de tentativas de melhoramento de uma solu��o
fes = POP_SIZE;               % Quantidade de avali��es da fun��o objetivo
%% pre alocacao da curva de converg�ncia
if SHOW_CONVERG
    %    MAX_ITER = floor((MAX_FES - POP_SIZE)/(2*POP_SIZE));
    %     fBestCurve = zeros(MAX_ITER +1,1);
    %     fesCurve = zeros(MAX_ITER +1,1);
    fBestCurve(1) = min(fobjValue);
    fesCurve(1) = fes;
end
iter = 1; % contador de itera��es
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
        
        % produzir nova solu��o
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
        
        % avaliar a nova solu��o
        fobjValueNew = fobj(xNew);
        fitValueNew = fitness(fobjValueNew);
        fes = fes +1;
        
        % quanto maior o fit melhor a solu��o
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
            
            % Gera nova solu��o
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
            % Avalia nova solu��o
            fobjValueNew = fobj(xNew);
            fitValueNew = fitness(fobjValueNew);
            fes = fes +1;
            
            % Quanto maior o fit melhor a solu��o
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
    % H� no m�ximo uma scout por fase
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