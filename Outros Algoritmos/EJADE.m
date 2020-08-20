function [xBest, fBest, fBestCurve, fesCurve] = EJADE(fobj, LB, UB, POP_SIZE, MAX_FES, seeConverg)
% Descrição
%     XXXX miniza a fobj usando a metaheurística XXXXX,
% conforme descrita em [1] e [2].
% Entradas:
%   fobj - Função objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada parâmetro
%   UB - Vetor linha com os limites superior de cada parâmetro
%   POP_SIZE - Inteiro com o tamanho da população
%   MAX_FES - Inteiro com o quantidade máxima de avalições da função objetivo
%   showConverg - Valor boleano que se for VERDADEIRO, ativará as saídas com os vetores 
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
%   [1] 
%   [2]
%% parâmetros do algoritmo
DIM = length(LB); % qtd de variaveis de design
c = 0.1; %*
p = 5;   %*

uF = 0.5;
uCR = 0.5;
POP_MIN = 4;
POP_MAX = POP_SIZE; %50
pop = POP_MAX;
%% Inicializa a populacao
xArchived = [];
x = LB + (UB - LB).*rand(pop, DIM);
fit = fobj(x);
fes = pop;
% if seeConverg
%     conver_curve = zeros(1, MAX_FES/pop);  população muda de tamanho,
%     cuidado!!!
% end
iter = 1;
if seeConverg
        fBestCurve(iter,1) = min(fit);
        fesCurve(iter,1) = fes;
end
while(fes + pop <= MAX_FES)
    S_CR = []; S_F = [];

    % Ajuste dos parametros
    CR = randNormal(uCR, 0.1, pop);% voltar aqui
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
        pbest  = id(randperm(ceil(p*pop/100),1));
        
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
    if length(xArchived) > POP_SIZE
        del = randperm(length(xArchived), length(xArchived) - POP_SIZE);
        xArchived(del,:) = [];
    end
    uCR = (1 - c)*uCR + c*mean(S_CR);
    uF = (1 - c)*uF + c*lehmarMean(S_F);
    fes = fes + pop; % number of objective function evaluation
    % redução da populacao
    newpop = floor((POP_MIN - POP_MAX)/MAX_FES*fes + POP_MAX);
    [~, id] = sort(fit,'descend');
    for i = 1:pop-newpop % pop reduction
        x(id(i),:) = [];
        fit(id(i)) = [];
    end
    pop = newpop;
    
    if seeConverg
        fBestCurve(iter,1) = min(fit);
        fesCurve(iter,1) = fes;
    end
    iter = iter+1;
end
[fBest, id] = min(fit);
xBest = x(id,:);
end

function xNew = boudaryCorrection(xNew, LB, UB, DIM, POP_SIZE)
%% LB e UB devem ser matrizes com dimensao [nBirds, dim]
u = (xNew < LB) | (xNew > UB);
randomMatrix = LB + (UB - LB).*rand(POP_SIZE, DIM);
xNew(u) = randomMatrix(u);
end

% function [x, fit, xArchived] =  updatePosition(x, fit, xNew, fitNew, xArchived, S_CR, CR, S_F, F)
% % isBetter = fitNew < fit;
% % if any(isBetter)
% %     fit(isBetter) = fitNew(isBetter);
% %     x(isBetter, :) = xNew(isBetter, :);
% %     idArq = ~isBetter;
% %     xArchived = [xArchived; x(idArq,:)];
% %
% % end
% if fitNew < fit
%     xArchived = [xArchived; x];
%     x = xNew;
%     fit = fitNew;
%     S_CR = [S_CR, CR];
%     S_F = [S_F, F];
% end
% end

function y = randCauchy(u, c, N)
% generate random numbers within Cauchy distribuition
% u - location parameter
% c - scale parameter
y = u + c*tan(pi*(rand(1,N) - 1/2));
end

function y = randNormal(mu, sd, N)
% generate random numbers within Normal distribuition
% u - mean
% c - standard deviation
y = randn(1,N)*sd + mu;
end

function y = lehmarMean(Sf)
y = sum(Sf.^2)/sum(Sf);
end
