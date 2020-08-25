function [xBest, fBest, fBestCurve, fesCurve] = JADE(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
% Descrição
%     JADE minimiza a fobj usando a metaheurística "Adaptive Differential 
% Evolution with Optional External Archive" conforme descrito em [1].
%
% Entradas:
%   fobj - Função objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada parâmetro
%   UB - Vetor linha com os limites superior de cada parâmetro
%   PARAM - Estrutura com o seguintes campos:
%      pop - Tamanho da população
%      p -  p% melhores, determines the greediness of the mutation strategy
%      c - controls the rate of parameter adaptation
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
%   [1] ZHANG, J.; SANDERSON, A. C. JADE: Adaptive differential evolution with optional external archive. IEEE Transactions on Evolutionary Computation, v. 13, n. 5, p. 945–958, 2009. 
%% parâmetros do algoritmo
% determines the greediness of the mutation strategy
p = PARAM.p; % os p% melhores

c = PARAM.c; % controls the rate of parameter adaptation
POP_SIZE = PARAM.pop; % tamanho da população
% valores iniciais para taxa de mutação e crossover
uF = 0.5;
uCR = 0.5;
%% Inicializa a populacao
DIM = length(LB); % qtd de variaveis de design
xArchived = [];
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
fit = fobj(x);
iter = 1;
fes = POP_SIZE;
if SHOW_CONVERG
    fBestCurve(iter) = min(fit);
    fesCurve(iter) = fes;
end
while(fes + POP_SIZE <= MAX_FES)
    S_CR = []; S_F = [];
    [~, id] = sort(fit);
    for i = 1:POP_SIZE
        %% 1 - Mutation (DE/current-to-pbest)
        % Ajuste dos parametros
        CR = randNormal(uCR, 0.1, 1);
        if CR > 1
            CR = 1;
        elseif CR < 0
            CR = 0;
        end
        
        F = randCauchy(uF, 0.1, 1);
        while F <= 0
            F = randCauchy(uF, 0.1, 1);
        end
        if F > 1
            F = 1;
        end

        % Selecionar um x dentre os p% melhores
        pbest  = id(randperm(round(p*POP_SIZE),1));
        
        % Selecionar x1 de P
        r1 = randperm(POP_SIZE, 2);
        r1 = r1(r1 ~= i);
        id1 = r1(1);
        x1 = x(id1,:);
        
        % Selecionar  x2 dentre geração atual e arquivada
        r2 = randperm(POP_SIZE + size(xArchived,1), 3);
        r2 = r2(r2 ~= i & r2 ~= id1);
        id2 = r2(1);
        
        if id2 <= POP_SIZE
            x2 = x(id2,:);
        else
            x2 = xArchived(id2-POP_SIZE,:);
        end
        
        % mutation vector v
        v(i,:) = x(i,:) + F.*(x(pbest,:) - x(i,:)) + F.*(x1 - x2);
        %% checa limites
        for d = 1:DIM
            if v(i,d) < LB(d)
                v(i,d) = (LB(d) + x(i,d))/2;
            elseif v(i,d) > UB(d)
                v(i,d) = (UB(d) + x(i,d))/2;
            end
        end
        %% 2 - Crossover
        jrand = randi(DIM);
        for j = 1:DIM
            if rand <= CR | j == jrand
                xNew(j) = v(i,j);
            else
                xNew(j) = x(i,j);
            end
        end
        %% 3 - Selection        
        % avalia a nova posicao
        fitNew = fobj(xNew); % avalia nova posicao
        if fitNew < fit(i)
            xArchived = [xArchived; x(i,:)];
            x(i,:) = xNew;
            fit(i) = fitNew;
            S_CR = [S_CR, CR];
            S_F = [S_F, F];
        end
    end
    fes = fes + POP_SIZE;
    % Randomly remove solutions from  xArchived
    if size(xArchived,1) > POP_SIZE
        del = randperm(size(xArchived,1), size(xArchived,1) - POP_SIZE);
        xArchived(del,:) = [];
    end
    uCR = (1 - c)*uCR + c*mean(S_CR);
    uF = (1 - c)*uF + c*lehmarMean(S_F);
    if SHOW_CONVERG
        fBestCurve(iter+1,1) = min(fit);
        fesCurve(iter+1,1) = fes;
    end
    iter = iter + 1;
end
[fBest, id] = min(fit);
xBest = x(id,:);
end

function y = randCauchy(u, c, N)
% Generate random numbers within Cauchy distribuition
% u - location parameter
% c - scale parameter
y = u + c*tan(pi*(rand(1,N) - 1/2));
end

function y = randNormal(mu, sd, N)
% Generate random numbers within Normal distribuition
% u - mean
% c - standard deviation
y = randn(1,N)*sd + mu;
end

function y = lehmarMean(Sf)
% Calcula a media de Lehmar
y = sum(Sf.^2)/sum(Sf);
end
