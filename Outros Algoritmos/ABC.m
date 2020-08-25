function [xBest, fBest, fBestCurve, fesCurve] = ABC(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
% Descrição
%   ABC minimiza a fobj usando a metaheurística "Artificial
%  Bee Colony", conforme descrito em [1] e [2]. O parâmetro limit é
%  especificado na ref [1]. O tratamento das restrições é adotado conforme
%  [2].
%
% Entradas:
%   fobj - Função objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada parâmetro
%   UB - Vetor linha com os limites superior de cada parâmetro
%   PARAM - Estrutura com o seguintes campos:
%      pop - Tamanho da população
%      limit - Numero de tentativas de melhoramento de uma fonte de alimentos
%   MAX_FES - Inteiro com o quantidade máxima de avalições da função objetivo
%   SHOW_CONVERG - Valor boleano que se for VERDADEIRO, ativará as saídas com os vetores 
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
%   [1] OLIVA, D.; CUEVAS, E.; PAJARES, G. Parameter identification of solar cells using artificial bee colony optimization. Energy, v. 72, p. 93 – 102, 2014.
%   [2] KARABOGA, D.; AKAY, B. A comparative study of Artificial Bee Colony algorithm. Applied Mathematics and Computation, Elsevier Inc., v. 214, n. 1, p. 108–132, 2009. ISSN 00963003. Disponível em: <http://dx.doi.org/10.1016/j.amc.2009.03.090>
%% parâmetros do algoritmo
DIM = length(LB); % qtd de variaveis de design

% numero de tentativas de melhoramento de uma fonte de alimentos
LIMIT = PARAM.limit;  % ref [1]
% A quantidade de fontes de alimento é igual ao quantidade de employed bees
NUM_FOODS = PARAM.pop;

% inicialização da população
x = LB + (UB - LB).*rand(NUM_FOODS, DIM); % População
fobjValue = fobj(x);           % Valor da função objetivo para cada solução candidata
fitValue = fitness(fobjValue); % Fitness de cada x. Quanto maior melhor
trial = zeros(1,NUM_FOODS);    % Quantidade de tentativas de melhoramento de uma solução
fes = NUM_FOODS;               % Quantidade de avalições da função objetivo
%% pre alocacao da curva de convergência
if SHOW_CONVERG
%    MAX_ITER = floor((MAX_FES - POP_SIZE)/(2*POP_SIZE));
%     fBestCurve = zeros(MAX_ITER +1,1);
%     fesCurve = zeros(MAX_ITER +1,1);
    fBestCurve(1) = min(fobjValue);
    fesCurve(1) = fes;
end
iter = 1; % contador de iterações
while(fes + 2*NUM_FOODS + 1 <= MAX_FES)
    %% Fase Employed bees
    for i = 1:NUM_FOODS
        % Parametro que será alterado
        j = randi(DIM);
        
        % selecionar uma abelha diferente de i
        k = randi(NUM_FOODS);   
        while(k == i)
            k = randi(NUM_FOODS);
        end
        
        % numero aleatorio no intevalo [-1, 1]
        phi = 2*rand - 1;   
        
        % produzir nova solução
        v = x(i,:);
        v(j) = x(i,j) + phi*(x(k,j) - x(i,j)); 
        
        % checar limitantes
        if v(j) < LB(j)
            v(j) = LB(j);
        elseif v(j) > UB(j)
            v(j) = UB(j);
        end
        
        % avaliar a nova solução
        fobjValueNew = fobj(v);
        fitValueNew = fitness(fobjValueNew);
        fes = fes +1;
        
        % quanto maior o fit melhor a solução
        if (fitValueNew > fitValue(i))
            x(i,:) = v;
            fitValue(i) = fitValueNew;
            fobjValue(i) = fobjValueNew;
            trial(i) = 0;
        else
            trial(i) = trial(i) + 1;
        end
    end

    %% Fase das onlooker bees
    prob = fitValue/sum(fitValue); % probabilidade de uma fonte de comida ser selecionada
    for i = 1:NUM_FOODS
        if(prob(i) > rand)
            j = randi(DIM);         % Parâmetro que será alterado
            k = randi(NUM_FOODS);   % Seleciona outra fonte de alimentos
            while(k == i)
                k = randi(NUM_FOODS);
            end
            
            phi = 2*rand - 1;
            
            % Gera nova solução
            v = x(i,:);
            v(j) = x(i,j) + phi*(x(k,j) - x(i,j));
            
            % Checa limitantes
            if v(j) < LB(j)
                v(j) = LB(j);
            elseif v(j) > UB(j)
                v(j) = UB(j);
            end
            
            % Avalia nova solução
            fobjValueNew = fobj(v);
            fitValueNew = fitness(fobjValueNew);
            fes = fes +1;
            

            % Quanto maior o fit melhor a solução
            if (fitValueNew > fitValue(i))
                x(i,:) = v;
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
        x(ind,:) = LB + (UB-LB).*rand(1,DIM); % ref [1]
        fobjValue(ind) = fobj(x(ind,:));
        fitValue(ind) = fitness(fobjValue(ind));
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
id = (fobjValue < 0);
fitval(id)  = 1 + abs(fobjValue(id));
end