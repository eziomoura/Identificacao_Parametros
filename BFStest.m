function [xBest, fBest, fBest_curve, fes_curve] = BFStest(fobj, LB, UB, POP_SIZE, MAX_FES, SHOW_CONVERG)
% Descrição
%     BFS minimiza a fobj usando a metaheurística "Birds Foraging Search"
% conforme descrito em [1].
%
% Entradas:
%   fobj - Função objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada parâmetro
%   UB - Vetor linha com os limites superior de cada parâmetro
%   POP_SIZE - Inteiro com o tamanho da população
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
%   [1] ZHANG, Z.; HUANG, C.; DONG, K.; HUANG, H. Birds foraging search: a novel population-based algorithm for global optimization. Memetic Computing, v. 11, n. 3, p. 221–250, 2019. 
%% parâmetros do algoritmo
nBirds = POP_SIZE; % tamanho da pop
%% Inicializa a populacao
DIM = length(LB); % qtd de variaveis de design
x = LB + (UB - LB).*rand(nBirds, DIM);
fit = fobj(x);    % fitness de cada individuo
fes = POP_SIZE;   % quantidade de avaliação da função objetivo
if SHOW_CONVERG
    MAX_ITER = floor((MAX_FES - POP_SIZE)/(3*POP_SIZE)); % numero maximo de iterações
    fBest_curve = zeros(MAX_ITER + 1, 1);
    fes_curve =  zeros(MAX_ITER + 1, 1);
    fBest_curve(1) = min(fit);
    fes_curve(1) = fes;
else 
    fBest_curve = NaN; fes_curve = NaN;
end
xIncNew = zeros(nBirds-1, DIM);
iter = 1;
POP_MIN = 6;
POP_MAX = POP_SIZE;
pop = POP_SIZE;
while(fes + 3*POP_SIZE <= MAX_FES)
    xOld = x;
    %% 1 - Flying search behavior phase
    [~, idBest] = min(fit);
    xBest = x(idBest,:);
    r1 = 2*rand(pop,1) - 1; % numeros aletorios entre -1 e 1
    xNew = abs(x - xBest).*exp(r1).*cos(2*pi*r1) + xBest; %equation (2)
    
    % verificar limites
    xNew = boudaryCorrection(xNew, LB, UB, DIM, pop);
    
    % avalia a nova posicao
    fitNew = fobj(xNew); % avalia nova posicao
    [x, fit] =  updatePosition(x, fit, xNew, fitNew);
    
    %% 2 - Execute territorial behavior phase
    %  * Territorial birds
    % identificar o territorial (melhor) os intrusos (subotimo, e os restantes)
    % O switch role mechanism é realizado aqui
    [~, idf] = sort(fit);
    xBest = x(idf(1),:);
    fitBest = fit(idf(1));
    
    xSecondBest = x(idf(2),:);
    xInc = x(idf(2:end),:);
    fitInc = fit(idf(2:end));
    
    r2 = 2*rand - 1; % numero aleatorio entre -1 e 1;
    lambda = xBest - xSecondBest; %lambda = 0.1*(UB(1,:) - LB(1,:));
    
    % gera novo xBest
    xBestNew = xBest + r2*lambda; % equation (3)
    
    % checa limites
    xBestNew = boudaryCorrection(xBestNew, LB, UB, DIM, 1);
    
    % avalia a nova posicao
    fitBestNew = fobj(xBestNew);
    [xBest, fitBest] =  updatePosition(xBest, fitBest, xBestNew, fitBestNew);
    
    % *incursion birds* xinc
    prob = 1 - iter/MAX_ITER; % Probabilidade de assustar os intrusos
    xIncNew = xInc;
    for i = 1:pop-1
        if rand >= prob
            IF = round(1 + rand);
            xIncNew(i,:) = xInc(i,:) + rand*(xBest - IF*xInc(i,:));
        else
            % generate random different integers;
            p = randperm(pop-1,5);
            p(p==i) = [];
            xIncNew(i,:) = xInc(i,:) + rand*(xInc(p(1),:) - xInc(p(2),:)) + rand*(xInc(p(3),:) - xInc(p(4),:));
        end
    end
    % Check boundary
    xIncNew = boudaryCorrection(xIncNew, LB, UB, DIM, pop-1);
    
    % avalia as novas posicões
    fitIncNew = fobj(xIncNew);
    [xInc, fitInc] = updatePosition(xInc, fitInc, xIncNew, fitIncNew);
    
    % Merge Population
%     x(idf(2:end),:) = xInc; fit(idf(2:end)) = fitInc;
%     x(idf(1),:) = xBest; fit(idf(1)) = fitBest;
    x(2:end,:) = xInc; fit(2:end) = fitInc;
    x(1,:) = xBest; fit(1) = fitBest;
%     
    %% Cognitive Behavior
    tol = 5*eps(x);
    isPosEqual = all(abs(x - xOld) < tol, 2);
    nEqual = sum(isPosEqual);
    zeta = (log(iter)/iter) * abs(xBest - rand(nEqual,1).*x(isPosEqual,:));
    xNew(isPosEqual,:) = randn(nEqual, 1).*zeta + xBest;
    xNew(~isPosEqual,:) = x(~isPosEqual,:) + rand(pop-nEqual, 1).*(x(~isPosEqual,:) - xOld(~isPosEqual,:));
    
    % Check boundary
    xNew = boudaryCorrection(xNew, LB, UB, DIM, pop);
    
    % avalia as novas posicões
    fitNew = fobj(xNew);
    [x, fit] = updatePosition(x, fit, xNew, fitNew);
    
    fes = fes + 3*pop;
    iter = iter + 1;
    
    % redução da populacao
    newpop = floor((POP_MIN - POP_MAX)/MAX_FES*fes + POP_MAX);
    [~, id] = sort(fit,'descend');
    
    % pop reduction
    x(id(1:pop-newpop),:) = [];
    fit(id(1:pop-newpop)) = [];
    pop = newpop;
    
    if SHOW_CONVERG
        fBest_curve(iter,1) = min(fit);
        fes_curve(iter,1) = fes;
    end
end
[fBest, id] = min(fit);
xBest = x(id,:);
end

function xNew = boudaryCorrection(xNew, LB, UB, DIM, nBirds)
% atribui valor aleatorio ao paramentro violar seus limites
    u = (xNew < LB) | (xNew > UB);
    randomMatrix = LB + (UB - LB).*rand(nBirds, DIM);
    xNew(u) = randomMatrix(u);
end

function [x, fit] =  updatePosition(x, fit, xNew, fitNew)
% realiza a greedy selection
    isBetter = fitNew < fit;
    if any(isBetter)
        fit(isBetter) = fitNew(isBetter);
        x(isBetter, :) = xNew(isBetter, :);
    end
end
