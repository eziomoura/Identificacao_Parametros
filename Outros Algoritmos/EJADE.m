function [xBest, fBest, fBestCurve, fesCurve] = EJADE(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
% Descrição
%     EJADE minimiza a fobj usando a metaheurística "Enhanced Adaptive 
% Differential Evolution" conforme descrita em [1].
% Autor não indicou como tratar as solution fora do range de busca. Aqui
% atribui-se o valor extremo mais próximo.
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
%   [1] LI, S.; GU, Q.; GONG, W.; NING, B. An enhanced adaptive differential evolution algorithm for parameter extraction of photovoltaic models. Energy Conversion and Management, v. 205, n. December 2019, p. 112443, 2020. 
%% parâmetros do algoritmo
c = PARAM.c; % os p% melhores
p = PARAM.p; % controls the rate of parameter adaptation
POP_MIN = PARAM.popMin;
POP_MAX = PARAM.popMax;
uF = 0.5;
uCR = 0.5;
%% Inicializa a populacao
pop = POP_MAX;
DIM = length(LB); % qtd de variaveis de design
xArchived = [];
x = LB + (UB - LB).*rand(pop, DIM);
fit = fobj(x);
fes = pop;
iter = 1;
if SHOW_CONVERG
    MAX_ITER = floor((MAX_FES - fes)/POP_MIN) + 1; % numero  maximo (com folga) de iterações
    fBestCurve = NaN(MAX_ITER, 1);
    fesCurve =  NaN(MAX_ITER, 1);
    fBestCurve(1) = min(fit);
    fesCurve(1) = fes;
else 
    fBest_curve = NaN; fes_curve = NaN;
end


while(fes + pop <= MAX_FES)
    S_CR = []; S_F = [];

    % Ajuste dos parametros
    CR = randNormal(uCR, 0.1, pop);
    F = randCauchy(uF, 0.1, pop);
    
    % Truncamentos
    for i = 1:pop
        if CR(i) > 1
            CR(i) = 1;
        elseif CR(i) < 0
            CR(i) = 0;
        end
        while F(i) <= 0
            F(i) = randCauchy(uF, 0.1, 1);
        end
        if F(i) > 1
            F(i) = 1;
        end
    end
    % Crossover rate sorting mechanism
     [~, id] = sort(fit); % eq 26
    CR = sort(CR); % eq 25
    CR(id) = CR;
    for i = 1:pop
        %% 1 - Mutation (DE/current-to-pbest)
        % Selecionar um x dentre os p% melhores
        pbest  = id(randperm(ceil(p*pop),1));
        
        % Selecionar x1 de P
        r1 = randperm(pop, 2);
        r1 = r1(r1 ~= i);
        id1 = r1(1);
        x1 = x(id1,:);
        
        % Selecionar  x2 de entre geração atual e arquivada
        r2 = randperm(pop + size(xArchived,1), 3);
        r2 = r2(r2 ~= i & r2 ~= id1);
        id2 = r2(1);
        if id2 <= pop
            x2 = x(id2,:);
        else
            x2 = xArchived(id2-pop,:);
        end
        % mutation vector v
        v(i,:) = x(i,:) + F(i).*(x(pbest,:) - x(i,:)) + F(i).*(x1 - x2);
        
        
        %% 2 - Crossover
        jrand = randi(DIM);
        for j = 1:DIM
            if rand <= CR(i) | j == jrand
                u(i,j) = v(i,j);
            else
                u(i,j) = x(i,j);
            end
        end
        %% 3 - Selection
        % verificar limites
        xNew = boudaryCorrection(u(i,:), LB, UB, DIM, 1);
        
        % avalia a nova posicao
        fitNew = fobj(xNew); % avalia nova posicao
        if fitNew < fit(i)
            xArchived = [xArchived; x(i,:)];
            x(i,:) = xNew;
            fit(i) = fitNew;
            S_CR = [S_CR, CR(i)];
            S_F = [S_F, F(i)];
        end
    end
    % Randomly remove solutions from  xArchived
    if size(xArchived,1) > POP_MAX
        del = randperm(size(xArchived,1), size(xArchived,1) - POP_MAX);
        xArchived(del,:) = [];
    end
    uCR = (1 - c)*uCR + c*mean(S_CR);
    uF = (1 - c)*uF + c*lehmarMean(S_F);
    fes = fes + pop; % number of objective function evaluation
    % redução da populacao
    newpop = floor(fes*(POP_MIN - POP_MAX)/MAX_FES + POP_MAX);
    [~, id] = sort(fit,'descend');
    
    % pop reduction
    x(id(1:pop-newpop),:) = [];
    fit(id(1:pop-newpop)) = [];
    pop = newpop;
    
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

% function xNew = boudaryCorrection(xNew, LB, UB, DIM, POP_SIZE)
% u = (xNew < LB) | (xNew > UB);
% randomMatrix = LB + (UB - LB).*rand(POP_SIZE, DIM);
% xNew(u) = randomMatrix(u);
% end

function xNew = boudaryCorrection(xNew, LB, UB, DIM, POP_SIZE)
u = xNew < LB;
xNew(u) = LB(u);
u = xNew > UB;
xNew(u) = UB(u);
end

function y = randCauchy(u, c, N)
% Generate random numbers within Cauchy distribuition
% u - location parameter
% c - scale parameter
y = u + c*tan(pi*(rand(1,N) - 1/2));
end

function y = randNormal(mu, sd, N)
% Generate random numbers within Normal distribuition
% mu - mean
% sd - standard deviation
y = randn(1,N)*sd + mu;
end

function y = lehmarMean(Sf)
% Media de Lehmar
y = sum(Sf.^2)/sum(Sf);
end
