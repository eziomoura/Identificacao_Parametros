function [xBest, fBest, fBestCurve, fesCurve] = JADE(fobj, LB, UB, POP_SIZE, MAX_FES, SHOW_CONVERG)
% Descri��o
%     XXXX miniza a fobj usando a metaheur�stica XXXXX,
% conforme descrita em [1] e [2].
% Entradas:
%   fobj - Fun��o objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada par�metro
%   UB - Vetor linha com os limites superior de cada par�metro
%   POP_SIZE - Inteiro com o tamanho da popula��o
%   MAX_FES - Inteiro com o quantidade m�xima de avali��es da fun��o objetivo
%   SHOW_CONVERG - Valor boleador que se for VERDADEIRO, ativar� as sa�das com os vetores 
%       referentes a curva de converg�ngia (converg_RMSE e converg_fes)
%        
% Sa�das:
%   xBest - Vetor com os par�metros que minimizam fobj
%   fBest - Valor da fobj avaliada em xBest
%   fBestCurve - Vetor com o fBest ao final de cada itera��o
%   fesCurve - Vetor com o n�mero de avali��es  da fun��o objetivo ao
%       final de cada itera��o
%
% Fontes:
%   [1] 
%   [2]
%% par�metros do algoritmo
DIM = length(LB); % qtd de variaveis de design
c = 0.1;
p = 5;

uF = 0.5;
uCR = 0.5;
%% Inicializa a populacao
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
        if CR>1
            CR = 1;
        elseif CR<0
            CR = 0;
        end
        
        F = randCauchy(uF, 0.1, 1);
        while F<=0
            F = randCauchy(uF, 0.1, 1);
        end
        if F > 1
            F = 1;
        end

        % Selecionar um x dentre os p% melhores
        pbest  = id(randperm(round(p*POP_SIZE/100),1));
        
        % Selecionar x1 de P
        r1 = randperm(POP_SIZE, 2);
        r1 = r1(r1 ~= i);
        id1 = r1(1);
        x1 = x(id1,:);
        
        % Selecionar  x2 dentre gera��o atual e arquivada
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
        
        %% 2 - Crossover
        jrand = randi(DIM);
        for j = 1:DIM
            if rand <= CR | j == jrand
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
            S_CR = [S_CR, CR];
            S_F = [S_F, F];
        end
    end
    fes = fes + POP_SIZE;
    % Randomly remove solutions from  xArchived
    if length(xArchived) > POP_SIZE
        del = randperm(length(xArchived), length(xArchived) - POP_SIZE);
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
