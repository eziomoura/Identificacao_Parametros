function [xBest, fBest, fBestCurve, fesCurve] = ABC(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
% Descri��o
%   ABC minimiza a fobj usando a metaheur�stica "Artificial
%  Bee Colony", conforme descrito em [1] e [2]. O par�metro limit �
%  especificado na ref [1]. O tratamento das restri��es � adotado conforme
%  [2].
%
% Entradas:
%   fobj - Fun��o objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada par�metro
%   UB - Vetor linha com os limites superior de cada par�metro
%   PARAM - Estrutura com o seguintes campos:
%      pop - Tamanho da popula��o
%      limit - Numero de tentativas de melhoramento de uma fonte de alimentos
%   MAX_FES - Inteiro com o quantidade m�xima de avali��es da fun��o objetivo
%   SHOW_CONVERG - Valor boleano que se for VERDADEIRO, ativar� as sa�das com os vetores 
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
%   [1] OLIVA, D.; CUEVAS, E.; PAJARES, G. Parameter identification of solar cells using artificial bee colony optimization. Energy, v. 72, p. 93 � 102, 2014.
%   [2] KARABOGA, D.; AKAY, B. A comparative study of Artificial Bee Colony algorithm. Applied Mathematics and Computation, Elsevier Inc., v. 214, n. 1, p. 108�132, 2009. ISSN 00963003. Dispon�vel em: <http://dx.doi.org/10.1016/j.amc.2009.03.090>
%% par�metros do algoritmo
DIM = length(LB); % qtd de variaveis de design

% numero de tentativas de melhoramento de uma fonte de alimentos
LIMIT = PARAM.limit;  % ref [1]
% A quantidade de fontes de alimento � igual ao quantidade de employed bees
NUM_FOODS = PARAM.pop;

% inicializa��o da popula��o
x = LB + (UB - LB).*rand(NUM_FOODS, DIM); % Popula��o
fobjValue = fobj(x);           % Valor da fun��o objetivo para cada solu��o candidata
fitValue = fitness(fobjValue); % Fitness de cada x. Quanto maior melhor
trial = zeros(1,NUM_FOODS);    % Quantidade de tentativas de melhoramento de uma solu��o
fes = NUM_FOODS;               % Quantidade de avali��es da fun��o objetivo
%% pre alocacao da curva de converg�ncia
if SHOW_CONVERG
%    MAX_ITER = floor((MAX_FES - POP_SIZE)/(2*POP_SIZE));
%     fBestCurve = zeros(MAX_ITER +1,1);
%     fesCurve = zeros(MAX_ITER +1,1);
    fBestCurve(1) = min(fobjValue);
    fesCurve(1) = fes;
end
iter = 1; % contador de itera��es
while(fes + 2*NUM_FOODS + 1 <= MAX_FES)
    %% Fase Employed bees
    for i = 1:NUM_FOODS
        % Parametro que ser� alterado
        j = randi(DIM);
        
        % selecionar uma abelha diferente de i
        k = randi(NUM_FOODS);   
        while(k == i)
            k = randi(NUM_FOODS);
        end
        
        % numero aleatorio no intevalo [-1, 1]
        phi = 2*rand - 1;   
        
        % produzir nova solu��o
        v = x(i,:);
        v(j) = x(i,j) + phi*(x(k,j) - x(i,j)); 
        
        % checar limitantes
        if v(j) < LB(j)
            v(j) = LB(j);
        elseif v(j) > UB(j)
            v(j) = UB(j);
        end
        
        % avaliar a nova solu��o
        fobjValueNew = fobj(v);
        fitValueNew = fitness(fobjValueNew);
        fes = fes +1;
        
        % quanto maior o fit melhor a solu��o
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
            j = randi(DIM);         % Par�metro que ser� alterado
            k = randi(NUM_FOODS);   % Seleciona outra fonte de alimentos
            while(k == i)
                k = randi(NUM_FOODS);
            end
            
            phi = 2*rand - 1;
            
            % Gera nova solu��o
            v = x(i,:);
            v(j) = x(i,j) + phi*(x(k,j) - x(i,j));
            
            % Checa limitantes
            if v(j) < LB(j)
                v(j) = LB(j);
            elseif v(j) > UB(j)
                v(j) = UB(j);
            end
            
            % Avalia nova solu��o
            fobjValueNew = fobj(v);
            fitValueNew = fitness(fobjValueNew);
            fes = fes +1;
            

            % Quanto maior o fit melhor a solu��o
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
    % H� no m�ximo uma scout por fase
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