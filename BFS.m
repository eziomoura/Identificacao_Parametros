function [xBest, fBest, fBest_curve, fes_curve] = BFS(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
% Descri��o
%     BFS minimiza a fobj usando a metaheur�stica "Birds Foraging Search"
% conforme descrito em [1].
%
% Entradas:
%   fobj - Fun��o objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada par�metro
%   UB - Vetor linha com os limites superior de cada par�metro
%   PARAM - Estrutura com o seguinte campo:
%      pop - Tamanho da popula��o
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
%   [1] ZHANG, Z.; HUANG, C.; DONG, K.; HUANG, H. Birds foraging search: a novel population-based algorithm for global optimization. Memetic Computing, v. 11, n. 3, p. 221�250, 2019. 
%% par�metros do algoritmo
POP_SIZE = PARAM.pop;
%% Inicializa a populacao
DIM = length(LB); % qtd de variaveis de design
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
fit = fobj(x);    % fitness de cada individuo
fes = POP_SIZE;   % quantidade de avalia��o da fun��o objetivo
if SHOW_CONVERG
    % max_iter  = #fes_paraUsarNoLoop/#quantidade usadada em cada itera��o
    % + a itera��o 0
    MAX_ITER = floor((MAX_FES - fes)/(3*POP_SIZE)) + 1; % numero maximo de itera��es
    fBest_curve = zeros(MAX_ITER, 1);
    fes_curve =  zeros(MAX_ITER, 1);
    fBest_curve(1) = min(fit);
    fes_curve(1) = fes;
else 
    fBest_curve = NaN; fes_curve = NaN;
end
xIncNew = zeros(POP_SIZE-1, DIM);
iter = 1;
while(fes + 3*POP_SIZE <= MAX_FES)
    xOld = x;
    %% 1 - Flying search behavior phase
    [~, idBest] = min(fit);
    xBest = x(idBest,:);
    r1 = 2*rand(POP_SIZE,1) - 1; % numeros aletorios entre -1 e 1
    xNew = abs(x - xBest).*exp(r1).*cos(2*pi*r1) + xBest; %equation (2)
    
    % verificar limites
    xNew = boudaryCorrection(xNew, LB, UB, DIM, POP_SIZE);
    
    % avalia a nova posicao
    fitNew = fobj(xNew); % avalia nova posicao
    [x, fit] =  updatePosition(x, fit, xNew, fitNew);
    
    %% 2 - Execute territorial behavior phase
    %  * Territorial birds
    % identificar o territorial (melhor) os intrusos (subotimo, e os restantes)
    % O switch role meachnism � realizado aqui
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
    for i = 1:POP_SIZE-1
        if rand >= prob
            IF = round(1 + rand);
            xIncNew(i,:) = xInc(i,:) + rand*(xBest - IF*xInc(i,:));
        else
            % generate random different integers;
            p = randperm(POP_SIZE-1,5);
            p(p==i) = [];
            xIncNew(i,:) = xInc(i,:) + rand*(xInc(p(1),:) - xInc(p(2),:)) + rand*(xInc(p(3),:) - xInc(p(4),:));
        end
    end
    % Check boundary
    xIncNew = boudaryCorrection(xIncNew, LB, UB, DIM, POP_SIZE-1);
    
    % avalia as novas posic�es
    fitIncNew = fobj(xIncNew);
    [xInc, fitInc] = updatePosition(xInc, fitInc, xIncNew, fitIncNew);
    
    % Merge Population
    x(idf(2:end),:) = xInc; fit(idf(2:end)) = fitInc;
    x(idf(1),:) = xBest; fit(idf(1)) = fitBest;
%     x(2:end,:) = xInc; fit(2:end) = fitInc;
%     x(1,:) = xBest; fit(1) = fitBest;
    
    %% Cognitive Behavior
    tol = 100*(eps(x) + eps(xOld));
    isPosEqual = all(abs(x - xOld) < tol, 2);
    nEqual = sum(isPosEqual);
    zeta = (log(iter)/iter) * abs(xBest - rand(nEqual,1).*x(isPosEqual,:));
    xNew(isPosEqual,:) = randn(nEqual, 1).*zeta + xBest;
    xNew(~isPosEqual,:) = x(~isPosEqual,:) + rand(POP_SIZE-nEqual, 1).*(x(~isPosEqual,:) - xOld(~isPosEqual,:));
    
    % Check boundary
    xNew = boudaryCorrection(xNew, LB, UB, DIM, POP_SIZE);
    
    % avalia as novas posic�es
    fitNew = fobj(xNew);
    [x, fit] = updatePosition(x, fit, xNew, fitNew);
    
    fes = fes + 3*POP_SIZE;
    iter = iter + 1;
    
    if SHOW_CONVERG
        fBest_curve(iter,1) = min(fit);
        fes_curve(iter,1) = fes;
    end
end
[fBest, id] = min(fit);
xBest = x(id,:);
end

function xNew = boudaryCorrection(xNew, LB, UB, DIM, nBirds)
% atribui valor aleatorio ao paramentro que violar seus limites
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
