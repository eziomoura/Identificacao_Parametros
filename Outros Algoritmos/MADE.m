function [xBest, fBest, fBestCurve, fesCurve] = MADE(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
% Descrição
%     MADE minimiza a fobj usando a metaheurística "memetic adaptive
%     differential evolution" conforme descrito em [1]. Em [2] fica
%     explícito que o autor utilizou o mesmo esquema de tratamento das
%     soluções que em SHADE
%
% Entradas:
%   fobj - Função objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada parâmetro
%   UB - Vetor linha com os limites superior de cada parâmetro
%   PARAM - Estrutura com o seguintes campos:
%      pop - Tamanho da população
%      p -  p% melhores, determines the greediness of the mutation strategy
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
% [1] LI, S.; GONG, W.; YAN, X.; HU, C.; BAI, D.; WANG, L. Parameter estimation of photovoltaic models with memetic adaptive differential evolution. Solar Energy, v. 190, n. September 2018, p. 465–474, 2019. 
% [2] https://github.com/wewnyin/wenyingong/blob/master/Software/MADE-code.zip
%%
DIM = length(LB); % dimensão do problema
%% parâmetros do algoritmo
% determines the greediness of the mutation strategy
H = PARAM.H;
POP_SIZE = PARAM.pop; % tamanho da população
epsilon =  PARAM.epsilon; % criterio para aplicar NM: f(x) < epsilon
NM_MAXFES =  200*DIM;% Máximo de avaliações da fobj dentro do NM
tol = 10^-8;
% valores iniciais para taxa de mutação e crossover
uF = repmat(0.5,H,1);
uCR = repmat(0.5,H,1);
%% Inicializa a populacao
xArchived = [];
fitArchived = [];
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
fit = fobj(x);
iter = 0; 
k = 1;
fes = POP_SIZE;

if SHOW_CONVERG
    MAX_ITER = floor((MAX_FES - fes)/POP_SIZE) + 1;
    fBestCurve = NaN(MAX_ITER,1);
    fesCurve = NaN(MAX_ITER,1);
    fBestCurve(1) = min(fit);
    fesCurve(1) = fes;
else
    fBestCurve = [];
    fesCurve = [];
end

while(fes + POP_SIZE <= MAX_FES)
    S_CR = []; S_F = []; df =[];
    [~, id] = sort(fit);
    for i = 1:POP_SIZE
        %% 1 - Mutation (DE/current-to-pbest)
        % Ajuste dos parametros
        r = randi(H);
        CR = randNormal(uCR(r), 0.1, 1);
        if CR > 1
            CR = 1;
        elseif CR < 0
            CR = 0;
        end
        
        F = randCauchy(uF(r), 0.1, 1);
        while F <= 0
            F = randCauchy(uF(r), 0.1, 1);
        end
        if F > 1
            F = 1;
        end

        % Selecionar um x dentre os p% melhores
        p = (0.2 - 2/POP_SIZE)*rand + 2/POP_SIZE;
        pbest  = id(randperm(round(p*POP_SIZE),1));
        
        
        % Selecionar x1 diferente i
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
        
        % checa limites de v
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
            df = [df, fit(i) - fitNew]; % array fitness value improvement
            xArchived = [xArchived; x(i,:)];
            fitArchived = [fitArchived; fit(i)];
            x(i,:) = xNew;
            fit(i) = fitNew;
            S_CR = [S_CR, CR];
            S_F = [S_F, F];
             
        end
    end
    fes = fes + POP_SIZE;
    %% Aplica Nelder Mead
    [~, idBest] = min(fit);
    if fit(idBest) < epsilon
        xBest = x(idBest,:);
        [x(idBest,:), fit(idBest), ~, fes_NM] = NelderMead(fobj, xBest, tol, NM_MAXFES, 1);
        fes = fes + fes_NM(end);
    end
    % Remove as piores solução de xArchived
    if size(xArchived,1) > POP_SIZE
        num2del = size(xArchived,1) - POP_SIZE;
        [~, del] = sort(fitArchived(1:POP_SIZE),'descend');
        xArchived(del(1:num2del),:) = [];
        fitArchived(del(1:num2del)) = [];
    end
    
    if ~isempty(S_CR) % or S_F
        w = df/sum(df);
        uCR(k) = sum(w.*S_CR);
        uF(k)  = sum(w.*S_F.^2)/sum(w.*S_F);
        k = k+1;
        if k > H
            k=1;
        end
    end
        iter = iter + 1;
    if SHOW_CONVERG
        fBestCurve(iter,1) = min(fit);
        fesCurve(iter,1) = fes;
    end

end
% Remove NaNs
fBestCurve = fBestCurve(1:iter, 1);
fesCurve = fesCurve(1:iter, 1);
%
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

function y = weightedLehmarMean(Sf)
% Calcula a media de Lehmar
y = sum(Sf.^2)/sum(Sf);
end


function y = lehmarMean(Sf)
% Calcula a media de Lehmar
y = sum(Sf.^2)/sum(Sf);
end
