function [xBest, fBest, fBestCurve, fesCurve] = DE(fobj, LB, UB, POP_SIZE, MAX_FES, SHOW_CONVERG)
% Descri��o
%     DE minimiza a fobj usando a metaheur�stica "Differential Evolution",
% conforme descrita em [1]. Como em [1] n�o � explicitado a forma de 
% tratamento das resti��es, aqui ser� atribuido um valor aleat�rio ao
% par�metro que ultrapasse seu limite superior ou inferior.
%
% Entradas:
%   fobj - Fun��o objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada par�metro
%   UB - Vetor linha com os limites superior de cada par�metro
%   POP_SIZE - Inteiro com o tamanho da popula��o
%   MAX_FES - Inteiro com o quantidade m�xima de avali��es da fun��o objetivo
%   SHOW_CONVERG - Valor boleano que se for VERDADEIRO, ativar� as sa�das com os vetores 
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
%   [1] STORN, R.; PRICE, K. Differential Evolution � A Simple and Efficient Heuristic for global Optimization over Continuous Spaces. Journal of Global Optimization, v. 11, n. 4, p. 341�359, 1 dez. 1997. 
%% par�metros do algoritmo
F = 0.95;
CR = 0.8;
%% Inicializa a populacao
DIM = length(LB); % qtd de variaveis de design
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
fit = fobj(x);
fes = POP_SIZE;
if SHOW_CONVERG
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
    
    if SHOW_CONVERG
        fBestCurve(iter) = min(fit);
        fesCurve(iter) = fes;
    end
end
[fBest, id] = min(fit);
xBest = x(id,:);
end

function xNew = boudaryCorrection(xNew, LB, UB, DIM, POP_SIZE)
% atribui valor aleatorio ao paramentro violar seus limites
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
