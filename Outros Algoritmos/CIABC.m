% OBS: Não obtive resultados significativante diferentes do ABC. 
% 1- há contradição entre o fluxograma e o texto, 
% 2- Autor não mostrou sua curva de convergia em relação a outros
% algoritmos
% 3- autor utilizou uma quantidade elevada de fes 

function [xBest, fBest, fBestCurve, fesCurve] = CIABC(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
% Descrição
%   CIABC minimiza a fobj usando a metaheurística "Chaotic Improved Artificial
%  Bee Colony", conforme descrito em [1]. Na fase 'Scout Bee' o autor diz
%  que as soluções que excederem o limite de tentativas de melhoramento
%  serão atualizadas usando a melhor solução ao invés de serem geradas
%  aleatoriamente. Porém em seu fluxograma é indicado o contrário. Aqui foi
%  usado a melhor solução.
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
%   [1] OLIVA, D.et al. A chaotic improved artificial bee colony for parameter estimation of photovoltaic cells. Energies, v. 10, n. 7, p. 1–19, 2017. ISSN 19961073
%% parâmetros do algoritmo
% numero de tentativas de melhoramento de uma fonte de alimentos
LIMIT = PARAM.limit;  % ref [1]
% A quantidade de fontes de alimento é igual ao quantidade de employed bees
POP_SIZE = PARAM.pop;

% inicialização da população
DIM = length(LB); % qtd de variaveis de design
x = LB + (UB - LB).*rand(POP_SIZE, DIM); % População inicial
fobjValue = fobj(x);           % Valor da função objetivo para cada solução candidata
fitValue = fitness(fobjValue); % Fitness de cada x. Quanto maior melhor
trial = zeros(1,POP_SIZE);    % Quantidade de tentativas de melhoramento de uma solução
fes = POP_SIZE;               % Quantidade de avalições da função objetivo


%% pre alocacao da curva de convergência
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
iter = 1; % contador de iterações
while(fes + 2*POP_SIZE <= MAX_FES)
    %% Fase Employed bees
    for i = 1:POP_SIZE
        % Parametro que será alterado
        j = randi(DIM);
        
        % selecionar uma abelha diferente de i
        k = randi(POP_SIZE);   
        while(k == i)
            k = randi(POP_SIZE);
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
    for i = 1:POP_SIZE
        if(prob(i) > CM(iter))
            j = randi(DIM);         % Parâmetro que será alterado
            k = randi(POP_SIZE);   % Seleciona outra fonte de alimentos
            while(k == i)
                k = randi(POP_SIZE);
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
    fes = fes + 2*POP_SIZE;
    %% SCOUT BEE
    % Há no máximo uma scout por fase
    % No CIABC é usada a melhor solução ao invés de uma solução aleatória

    [~, ind] =  max(trial);
    if (trial(ind) > LIMIT)
        [~, idBest] = min(fobjValue);
        trial(ind) = 0;
        x(ind,:) = x(idBest,:); % gera uma nova solução em torno da melhor.
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
