%% ver como fazer o clamp e o constrain position
function [xBest, fBest, fBestCurve, fesCurve] = ELPSO(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
% Descrição
%     ELPSO minimiza a fobj usando a metaheurística "Enhanced leader 
% particle swarm optimisation", conforme descrita em [1] e [2].
% Entradas:
%   fobj - Função objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada parâmetro
%   UB - Vetor linha com os limites superior de cada parâmetro
%   PARAM - Inteiro com o tamanho da população
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
%   [1] JORDEHI, A. R. Enhanced leader particle swarm optimisation (ELPSO): An efficient algorithm for parameter estimation of photovoltaic (PV) cells and modules. Solar Energy, v. 159, n. March 2017, p. 78–87, 2018.
%   [2] JORDEHI, A. R. Enhanced leader PSO (ELPSO): A new PSO variant for solving global optimisation problems. Applied Soft Computing Journal, v. 26, p. 401–417, 2015. 
%% parâmetros do algoritmo
c1 = PARAM.c1; % personal acceleration coefficient
c2 = PARAM.c2; % social acceleration coefficien
w = PARAM.w;   % inertia weight
h = PARAM.h;   % standard deviation (Gaussian)
s = PARAM.s;   % scale factor (Cauchy)
F = PARAM.F;   % scale factor of DE-based mutation
POP_SIZE = PARAM.pop; % population
%% inicialize positions and velocity
vMax = (UB - LB);
DIM = length(LB);
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
v = zeros(size(x));
fobjValue = fobj(x);
pBest = x; % personal best of each particle
pBestFit = fobjValue;
[gBestFit, idBest] = min(fobjValue);
gBest = x(idBest,:); % global best

%% pre alocacao
fes = POP_SIZE;
if SHOW_CONVERG
    MAX_ITER = floor((MAX_FES - POP_SIZE)/POP_SIZE);
    fBestCurve = zeros(MAX_ITER + 1, 1);
    fesCurve = zeros(MAX_ITER + 1, 1);
    fBestCurve(1) = gBestFit;
    fesCurve(1) = fes;
end
iter = 1;
while(fes+POP_SIZE <= MAX_FES)
        % calc inercia weight
        w = 0.9 -(0.5)/(iter/MAX_ITER);
        
    for i = 1:POP_SIZE
        % update velocity
        v(i,:) = w*v(i,:) + c1*rand*(pBest(i,:) - x(i,:))...
                          + c2*rand*(gBest - x(i,:));
        
        % clamp velocity
        v(i,:) = min(v(i,:),  vMax);
        v(i,:) = max(v(i,:), -vMax);
        
        % bound position
        x(i,:) = x(i,:) + v(i,:);
        x(i,:) = min(x(i,:), UB);
        x(i,:) = max(x(i,:), LB);
        

        
        % update personal best
        fobjValue(i) = fobj(x(i,:));
        if fobjValue(i) < pBestFit(i)
            pBest(i,:) = x(i,:);
            pBestFit(i) = fobjValue(i);
            
            % update global best (so far)
            if fobjValue(i) < gBestFit
                gBest = x(i,:);
                gBestFit = fobjValue(i);
            end
        end
    end
    fes = fes + POP_SIZE;
    
    % 1- Gaussian mutation
    gBestNew = gBest + (UB - LB) .* randNormal(0, h, DIM);
    gBestNewfit = fobj(gBestNew);
    if gBestNewfit < gBestFit
        gBest = gBestNew;
        gBestFit = gBestNewfit;
    end
    % update standard deviation
    h = h - 1/MAX_ITER;
    
    % 2- Cauchy mutation
    gBestNew = gBest + (UB - LB) .* randCauchy(0, s, DIM);
    gBestNewfit = fobj(gBestNew);
    if gBestNewfit < gBestFit
        gBest = gBestNew;
        gBestFit = gBestNewfit;
    end
    % update scale factor
    s = s - 1/MAX_ITER;
    
    % 3- Opposition-based mutation (seradamente aplicado a cada dim)
    for d=1:DIM
        gBestNew = gBest;
        gBestNew(d) = LB(d) + UB(d) - gBest(d);
        gBestNewfit = fobj(gBestNew);
        if gBestNewfit < gBestFit
            gBest = gBestNew;
            gBestFit = gBestNewfit;
        end
    end

    
    % 4- Opposition-based mutation (aplicado a todo o vetor)
    gBestNew = LB + UB - gBest;
    gBestNewfit = fobj(gBestNew);
    if gBestNewfit < gBestFit
        gBest = gBestNew;
        gBestFit = gBestNewfit;
    end
    
    % 5- DE-based mutation
    id = randperm(POP_SIZE,2);
    gBestNew = gBest + F*(pBest(id(1),:) - pBest(id(2),:));
    gBestNewfit = fobj(gBestNew);
    if gBestNewfit < gBestFit
        gBest = gBestNew;
        gBestFit = gBestNewfit;
    end
    fes = fes + DIM + 4;
    
    iter = iter + 1;
    if SHOW_CONVERG
        fBestCurve(iter,1) = gBestFit;
        fesCurve(iter,1) = fes;
    end
end
xBest = gBest;
fBest = gBestFit;
end

function xNew = boudaryCorrection(xNew, LB, UB, DIM, nBirds)
% atribui valor aleatorio ao paramentro violar seus limites
    u = (xNew < LB) | (xNew > UB);
    randomMatrix = LB + (UB - LB).*rand(nBirds, DIM);
    xNew(u) = randomMatrix(u);
end

function y = randNormal(mu, sd, N)
% Generate random numbers within Normal distribuition
% mu - mean
% sd - standard deviation
y = randn(1,N)*sd + mu;
end

function y = randCauchy(u, c, N)
% Generate random numbers within Cauchy distribuition
% u - location parameter
% c - scale parameter
y = u + c*tan(pi*(rand(1,N) - 1/2));
end