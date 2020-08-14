function [xBest, fBest, fBest_curve, fes_curve] = BFS(fobj, LB, UB, POP_SIZE, MAX_FES, SHOW_CONVERG)
% Descrição
%     XXXX miniza a fobj usando a metaheurística XXXXX,
% conforme descrita em [1] e [2].
% Entradas:
%   fobj - Função objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada parâmetro
%   UB - Vetor linha com os limites superior de cada parâmetro
%   POP_SIZE - Inteiro com o tamanho da população
%   MAX_FES - Inteiro com o quantidade máxima de avalições da função objetivo
%   showConverg - Valor boleador que se for VERDADEIRO, ativará as saídas com os vetores 
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
nBirds = POP_SIZE; % tamanho da pop
DIM = length(LB); % qtd de variaveis de design
%% Inicializa a populacao
x = LB + (UB - LB).*rand(nBirds, DIM);
fit = fobj(x);
fes = POP_SIZE; % quantidade de avaliação da função objetivo
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
while(fes + 3*POP_SIZE <= MAX_FES)
    xOld = x;
    %% 1 - Flying search behavior phase
    [~, idBest] = min(fit);
    xBest = x(idBest,:);
    r1 = 2*rand(nBirds,1) - 1; % numeros aletorios entre -1 e 1
    xNew = abs(x - xBest).*exp(r1).*cos(2*pi*r1) + xBest; %equation (2)
    
    % verificar limites
    xNew = boudaryCorrection(xNew, LB, UB, DIM, nBirds);
    
    % avalia a nova posicao
    fitNew = fobj(xNew); % avalia nova posicao
    [x, fit] =  updatePosition(x, fit, xNew, fitNew);
    
    %% 2 - Execute territorial behavior phase
    %  * Territorial birds
    % identificar o territorial (melhor) os intrusos (subotimo, e os restantes)
    % O switch role meachnism é realizado aqui
    [~, idf] = sort(fit);
    xBest = x(idf(1),:);
    fitBest = fit(idf(1));
    
    xSecondBest = x(idf(2),:);
    xInc = x(idf(2:end),:);
    fitInc = fit(idf(2:end));
    
    r2 = 2*rand - 1; % numero aleatorio entre -1 e 1;
    lambda = xBest - xSecondBest; %lambda = 0.1*(UB(1,:) - LB(1,:));
    xBestNew = xBest + r2.*lambda; % equation (3)
    xBestNew = boudaryCorrection(xBestNew, LB(1,:), UB(1,:), DIM, 1);
    fitBestNew = fobj(xBestNew);
    [xBest, fitBest] =  updatePosition(xBest, fitBest, xBestNew, fitBestNew);
    
    % *incursion birds* xinc
    prob = 1 - iter/MAX_ITER; % Probabilidade de assustar os intrusos
    for i = 1:nBirds-1
        if prob <= rand
            IF = round(1 + rand);
            xIncNew(i,:) = xInc(i,:) + rand*(xBest - IF*xInc(i,:));
        else
            % generate random different integers;
            p = randperm(nBirds-1,5);
            p(p==i) = [];
            %             temp = [1:(i-1), (i+1):(nBirds-1)];
            %             p = temp(randi(nBirds-2,[1, 4]));
%             d = randi(DIM);
            xIncNew(i,:) = xInc(i,:) + rand*(xInc(p(1),:) - xInc(p(2),:)) + rand*(xInc(p(3),:) - xInc(p(4),:));
        end
    end
    % Check boundary
    xIncNew = boudaryCorrection(xIncNew, LB, UB, DIM, nBirds-1);
    
    % avalia as novas posicões
    fitIncNew = fobj(xIncNew);
    [xInc, fitInc] = updatePosition(xInc, fitInc, xIncNew, fitIncNew);
    
    % Merge Population
%     x(idf(2:end),:) = xInc; fit(idf(2:end)) = fitInc;
%     x(idf(1),:) = xBest; fit(idf(1)) = fitBest;
    x(2:end,:) = xInc; fit(2:end) = fitInc;
    x(1,:) = xBest; fit(1) = fitBest;
    
    %% Cognitive Behavior
    isPosEqual = all(x == xOld, 2);
    nEqual = sum(isPosEqual);
    zeta = log(iter)/iter * abs(xBest - rand(nEqual,DIM).*x(isPosEqual,:));
    xNew(isPosEqual,:) = randn(nEqual, DIM).*zeta + xBest;
    xNew(~isPosEqual,:) = x(~isPosEqual,:) + rand(nBirds-nEqual, DIM).*(x(~isPosEqual,:) - xOld(~isPosEqual,:));
    
    % Check boundary
    xNew = boudaryCorrection(xNew, LB, UB, DIM, nBirds);
    
    % avalia as novas posicões
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
%% LB e UB devem ser matrizes com dimensao [nBirds, dim]
    u = (xNew < LB) | (xNew > UB);
    randomMatrix = LB + (UB - LB).*rand(nBirds, DIM);
    xNew(u) = randomMatrix(u);
end

function [x, fit] =  updatePosition(x, fit, xNew, fitNew)
    isBetter = fitNew < fit;
    if any(isBetter)
        fit(isBetter) = fitNew(isBetter);
        x(isBetter, :) = xNew(isBetter, :);
    end
end
