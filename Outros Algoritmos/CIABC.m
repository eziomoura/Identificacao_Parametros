% OBS: N�o obtive resultados significativante diferentes do ABC. 
% 1- h� contradi��o entre o fluxograma e o texto, 
% 2- Autor n�o mostrou sua curva de convergia em rela��o a outros
% algoritmos
% 3- autor utilizou uma quantidade elevada de fes 

function [xBest, fBest, fBestCurve, fesCurve] = CIABC(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
% Descri��o
%   CIABC minimiza a fobj usando a metaheur�stica "Chaotic Improved Artificial
%  Bee Colony", conforme descrito em [1]. Na fase 'Scout Bee' o autor diz
%  que as solu��es que excederem o limite de tentativas de melhoramento
%  ser�o atualizadas usando a melhor solu��o ao inv�s de serem geradas
%  aleatoriamente. Por�m em seu fluxograma � indicado o contr�rio. Aqui foi
%  usado a melhor solu��o.
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
%   [1] OLIVA, D.et al. A chaotic improved artificial bee colony for parameter estimation of photovoltaic cells. Energies, v. 10, n. 7, p. 1�19, 2017. ISSN 19961073
%% par�metros do algoritmo
% numero de tentativas de melhoramento de uma fonte de alimentos
LIMIT = PARAM.limit;  % ref [1]
% A quantidade de fontes de alimento � igual ao quantidade de employed bees
POP_SIZE = PARAM.pop;

% inicializa��o da popula��o
DIM = length(LB); % qtd de variaveis de design
x = LB + (UB - LB).*rand(POP_SIZE, DIM); % Popula��o inicial
fobjValue = fobj(x);           % Valor da fun��o objetivo para cada solu��o candidata
fitValue = fitness(fobjValue); % Fitness de cada x. Quanto maior melhor
trial = zeros(1,POP_SIZE);    % Quantidade de tentativas de melhoramento de uma solu��o
fes = POP_SIZE;               % Quantidade de avali��es da fun��o objetivo


%% pre alocacao da curva de converg�ncia
if SHOW_CONVERG
    MAX_ITER = floor((MAX_FES - POP_SIZE)/(2*POP_SIZE));
    fBestCurve = zeros(MAX_ITER +1,1);
    fesCurve = zeros(MAX_ITER +1,1);
    fBestCurve(1) = min(fobjValue);
    fesCurve(1) = fes;
end
%% mapa caotico de Tent
CM = tentMap(rand, MAX_ITER);
%%
iter = 1; % contador de itera��es
while(fes + 2*POP_SIZE <= MAX_FES)
    %% Fase Employed bees
    for i = 1:POP_SIZE
        % Parametro que ser� alterado
        j = randi(DIM);
        
        % selecionar uma abelha diferente de i
        k = randi(POP_SIZE);   
        while(k == i)
            k = randi(POP_SIZE);
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
    for i = 1:POP_SIZE
        if(prob(i) > CM(iter))
            j = randi(DIM);         % Par�metro que ser� alterado
            k = randi(POP_SIZE);   % Seleciona outra fonte de alimentos
            while(k == i)
                k = randi(POP_SIZE);
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
    fes = fes + 2*POP_SIZE;
    %% SCOUT BEE
    % H� no m�ximo uma scout por fase
    % No CIABC � usada a melhor solu��o ao inv�s de uma solu��o aleat�ria

    [~, ind] =  max(trial);
    if (trial(ind) > LIMIT)
        [~, idBest] = min(fobjValue);
        trial(ind) = 0;
        x(ind,:) = x(idBest,:); % gera uma nova solu��o em torno da melhor.
        fobjValue(ind) = fobjValue(idBest);
        fitValue(ind) = fitValue(idBest);
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
% calcula o fitness
fitval = zeros(size(fobjValue));
id = (fobjValue >= 0);
fitval(id) = 1./(fobjValue(id)+1);
id = (fobjValue < 0);
fitval(id)  = 1 + abs(fobjValue(id));
end

function g = tentMap(p, N)
% gera mapa caotico de tent com N valores e valor incial p
g = zeros(1,N);
g(1) = p;
for i = 2:N
    if g(i-1) < 0.7
        g(i) = g(i-1)/0.7;
    else
        g(i) = g(i-1)*(1 - g(i-1))/0.3;
    end
end
end
