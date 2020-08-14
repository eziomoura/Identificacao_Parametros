function [Iph, I0, n, Rs, Rp, RMSE, converg_RMSE, converg_fes] = ABC(fobj, LB, UB, POP_SIZE, MAX_FES, showConverg)
% Descri��o
%   ABC estima os par�metros el�tricos usando a metaheur�stica Artificial
% Bee Colony, conforme descrita em [1] e [2]. Nenhum dos autores([1] e [2])
% indicou sua abordagem em rela��o as restri��es do espa�o de busca.
% Aqui utilizamos a est�tegia de selecionar um valor alet�rio dentro do
% espa�o de busca para o par�metro que exceder seu limite inferior ou
% superior.
%
% Entradas:
%   fobj - a fun��o objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada par�metro
%   UB - Vetor linha com os limites superior de cada par�metro
%   POP_SIZE - Inteiro com o tamanho da popula��o
%   MAX_FES - Inteiro com o quantidade m�xima de avali��es da fun��o objetivo
%   showConverg - Valor boleador que se for VERDADEIRO, ativar� as sa�das com os vetores 
%       referentes a curva de converg�ngia (converg_RMSE e converg_fes)
%        
% Sa�das:
%   Iph - Corrente fotogerada [A]
%   I0  - Corrente de satura��o do diodo [A]
%   n   - �ndice de idealidade do diodo
%   Rs  - Resist�ncia S�rie [Ohms]
%   Rp  - Resist�ncia Paralelo [Ohms]
%   RMSE - Erro quadr�tico m�dio [A]
%   converg_RMSE - Vetor com o RMSE ao final de cada itera��o
%   converg_fes - Vetor com o n�mero de avali��es  da fun��o objetivo ao
%       final de cada itera��o
%
% Fontes:
%   [1] OLIVA, D.; CUEVAS, E.; PAJARES, G. Parameter identification of solar cells using artificial bee colony optimization. Energy, v. 72, p. 93 � 102, 2014.
%   [2] KARABOGA, D.; BASTURK, B. A powerful and efficient algorithm for numerical function optimization: Artificial bee colony (ABC) algorithm. Journal of Global Optimization, v. 39, n. 3,p. 459�471, 2007. ISSN 09255001
%% par�metros do algoritmo
DIM = length(LB); % qtd de variaveis de design
LIMIT = 100; % numero de tentativas de melhoramento de uma fonte de alimentos
% A quantidade de fontes de alimento � igual ao quantidade de employed bees
NUM_FOODS = POP_SIZE; % number of employed bees
%% inicializa��o da popula��o
x = LB + (UB - LB).*rand(NUM_FOODS, DIM); % Popula��o
fobjValue = fobj(x); % Valor da fun��o objetivo para cada solu��o candidata
fitValue = fitness(fobjValue); % Fitness de cada x. Quanto maior melhor.
trial = zeros(1,NUM_FOODS); % Quantidade de tentativas de melhoramento de uma solu��o
fes = NUM_FOODS; % Quantidade de avali��es da fun��o objetivo
%% pre alocacao da curva de converg�ncia
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
        % Parametro que ser� alterado
        j = randi(DIM);
        
        % seleccionar uma abelha diferente de i
        k = randi(NUM_FOODS);   
        while(k == i)
            k = randi(NUM_FOODS);
        end
        
        phi = 2*rand - 1; % random number between -1 and 1
        
        % produzir nova solu��o
        v = x(i,:);
        v(j) = x(i,j) + phi*(x(k,j) - x(i,j)); 
        
        % checar limitantes
        v(v(j) < LB(j) | v(j) > UB(j)) = UB(j) - rand*(UB(j) - LB(j));
        
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
            v(v(j) < LB(j) | v(j) > UB(j)) = UB(j) - rand*(UB(j) - LB(j));
            
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
    
    %% SCOUT BEE
    % H� no m�ximo uma scout por fase
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