function [Iph, I0, n, Rs, Rp, RMSE, converg_RMSE, converg_fes] = ABC(fobj, LB, UB, POP_SIZE, MAX_FES, showConverg)
% Descrição
%   ABC estima os parâmetros elétricos usando a metaheurística Artificial
% Bee Colony, conforme descrita em [1] e [2]. Nenhum dos autores([1] e [2])
% indicou sua abordagem em relação as restrições do espaço de busca.
% Aqui utilizamos a estátegia de selecionar um valor aletório dentro do
% espaço de busca para o parâmetro que exceder seu limite inferior ou
% superior.
%
% Entradas:
%   fobj - a função objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada parâmetro
%   UB - Vetor linha com os limites superior de cada parâmetro
%   POP_SIZE - Inteiro com o tamanho da população
%   MAX_FES - Inteiro com o quantidade máxima de avalições da função objetivo
%   showConverg - Valor boleador que se for VERDADEIRO, ativará as saídas com os vetores 
%       referentes a curva de convergêngia (converg_RMSE e converg_fes)
%        
% Saídas:
%   Iph - Corrente fotogerada [A]
%   I0  - Corrente de saturação do diodo [A]
%   n   - Índice de idealidade do diodo
%   Rs  - Resistência Série [Ohms]
%   Rp  - Resistência Paralelo [Ohms]
%   RMSE - Erro quadrático médio [A]
%   converg_RMSE - Vetor com o RMSE ao final de cada iteração
%   converg_fes - Vetor com o número de avalições  da função objetivo ao
%       final de cada iteração
%
% Fontes:
%   [1] OLIVA, D.; CUEVAS, E.; PAJARES, G. Parameter identification of solar cells using artificial bee colony optimization. Energy, v. 72, p. 93 – 102, 2014.
%   [2] KARABOGA, D.; BASTURK, B. A powerful and efficient algorithm for numerical function optimization: Artificial bee colony (ABC) algorithm. Journal of Global Optimization, v. 39, n. 3,p. 459–471, 2007. ISSN 09255001
%% parâmetros do algoritmo
DIM = length(LB); % qtd de variaveis de design
LIMIT = 100; % numero de tentativas de melhoramento de uma fonte de alimentos
% A quantidade de fontes de alimento é igual ao quantidade de employed bees
NUM_FOODS = POP_SIZE; % number of employed bees
%% inicialização da população
x = LB + (UB - LB).*rand(NUM_FOODS, DIM); % População
fobjValue = fobj(x); % Valor da função objetivo para cada solução candidata
fitValue = fitness(fobjValue); % Fitness de cada x. Quanto maior melhor.
trial = zeros(1,NUM_FOODS); % Quantidade de tentativas de melhoramento de uma solução
fes = NUM_FOODS; % Quantidade de avalições da função objetivo
%% pre alocacao da curva de convergência
if showConverg
    MAX_ITER = floor((MAX_FES - POP_SIZE)/POP_SIZE);
    converg_RMSE = zeros(MAX_ITER +1,1);
    converg_fes = zeros(MAX_ITER +1,1);
    converg_RMSE(1) = sqrt(min(fobjValue));
    converg_fes(1) = fes;
end
iter = 1;
while(fes + 2*NUM_FOODS + 1 <= MAX_FES)
    %% Fase Employed bees
    for i = 1:NUM_FOODS
        % Parametro que será alterado
        j = randi(DIM);
        
        % seleccionar uma abelha diferente de i
        k = randi(NUM_FOODS);   
        while(k == i)
            k = randi(NUM_FOODS);
        end
        
        phi = 2*rand - 1; % random number between -1 and 1
        
        % produzir nova solução
        v = x(i,:);
        v(j) = x(i,j) + phi*(x(k,j) - x(i,j)); 
        
        % checar limitantes
        v(v(j) < LB(j) | v(j) > UB(j)) = UB(j) - rand*(UB(j) - LB(j));
        
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
            v(v(j) < LB(j) | v(j) > UB(j)) = UB(j) - rand*(UB(j) - LB(j));
            
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
    
    %% SCOUT BEE
    % Há no máximo uma scout por fase
    fes = fes + 2*NUM_FOODS;
    [~, ind] =  max(trial);
    if (trial(ind) > LIMIT)
        trial(ind) = 0;
        x(ind,:) = LB + (UB-LB).*rand(1,DIM);
        fobjValue(ind) = fobj(x(ind,:));
        fitValue(ind) = fitness(fobjValue(ind));
        fes = fes +1;
    end
    % Atualizar curva de convergencia
    iter = iter +1;
    if showConverg
        converg_RMSE(iter,1) = sqrt(min(fobjValue));
        converg_fes(iter,1) = fes;
    end
end
[~, id] = max(fitValue);
bestGlobalMSE = fobjValue(id);
RMSE = sqrt(bestGlobalMSE);
gbest = x(id,:);
Iph = gbest(1);
I0 = gbest(2);
n = gbest(3);
Rs = gbest(4);
Rp = gbest(5);
end

function fitval = fitness(fobjValue)
fitval = zeros(size(fobjValue));
id = (fobjValue >= 0);
fitval(id) = 1./(fobjValue(id)+1);
id = (fobjValue < 0);
fitval(id)  = 1 + abs(fobjValue(id));
end