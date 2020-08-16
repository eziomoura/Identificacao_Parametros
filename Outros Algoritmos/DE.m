function [xBest, fBest, fBestCurve, fesCurve] = DE(fobj, LB, UB, POP_SIZE, MAX_FES, SHOW_CONVERG)
% Descrição
%     XXXX miniza a fobj usando a metaheurística XXXXX,
% conforme descrita em [1] e [2].
% Entradas:
%   fobj - Função objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada parâmetro
%   UB - Vetor linha com os limites superior de cada parâmetro
%   POP_SIZE - Inteiro com o tamanho da população
%   MAX_FES - Inteiro com o quantidade máxima de avalições da função objetivo
%   SHOW_CONVERG - Valor boleador que se for VERDADEIRO, ativará as saídas com os vetores 
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
F = 0.95;
CR = 0.8;
%maxIter = max_fes/(N); % numero maximo de iterações
%% Inicializa a populacao
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
fit = fobj(x);
fes = POP_SIZE;
if seeConverg
    MAX_ITER = (MAX_FES - POP_SIZE)/POP_SIZE;
    fBestCurve = zeros(MAX_ITER,1);
    fesCurve = zeros(MAX_ITER,1);
    fBestCurve(1) = min(fit);
    fesCurve(1) = fes;
end
iter = 1;
while(fes + POP_SIZE <= MAX_FES)
    %% 1 - Mutation (DE/RAND/1)
    for i = 1:POP_SIZE
        r = randperm(POP_SIZE, 3);
        r = r(r ~= i);
        v(i,:) = x(i,:) + F.*(x(r(1),:) - x(r(2),:));
    end
    %% 2 - Crossover
    for i = 1:POP_SIZE
        p = randi(DIM);
        for j = 1:DIM
            if rand <= CR | j == p
                u(i,j) = v(i,j);
            else
                u(i,j) = x(i,j);
            end
        end
    end
    %% 
    % verificar limites
    xNew = u;
    xNew = boudaryCorrection(xNew, LB, UB, DIM, POP_SIZE);
    
    % avalia a nova posicao
    fitNew = fobj(xNew); % avalia nova posicao
    [x, fit] =  updatePosition(x, fit, xNew, fitNew);
    
    fes = fes + POP_SIZE;
    iter = iter +1;
    
    if seeConverg
        fBestCurve(iter) = min(fit);
        fesCurve(iter) = fes;
    end
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

function [x, fit] =  updatePosition(x, fit, xNew, fitNew)
    isBetter = fitNew < fit;
    if any(isBetter)
        fit(isBetter) = fitNew(isBetter);
        x(isBetter, :) = xNew(isBetter, :);
    end
end
