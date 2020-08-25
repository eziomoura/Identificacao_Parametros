function [xBest, fBest, fBestCurve, fesCurve] = PSO(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
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
c1 = PARAM.c1; % personal acceleration coefficient
c2 = PARAM.c2; % social acceleration coefficien
w = PARAM.w;
vMax = (UB - LB);
POP_SIZE = PARAM.pop;
%% inicialize positions and velocity
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
    fBestCurve(1) = min(fobjValue);
    fesCurve(1) = fes;
end
iter = 1;
while(fes+POP_SIZE <= MAX_FES)
    for i = 1:POP_SIZE
        % update velocity
        v(i,:) = w*v(i,:) + c1*rand(1,DIM).*(pBest(i,:) - x(i,:)) + c2*rand(1,DIM).*(gBest - x(i,:));
        v(i,:) = min(v(i,:),  vMax);
        v(i,:) = max(v(i,:), -vMax);
        
        % update position
        x(i,:) = x(i,:) + v(i,:);
        x(i,:) = min(x(i,:), UB);
        x(i,:) = max(x(i,:), LB);
        
        % calc inercia weight
        w = 0.9 -(0.4)/(iter/MAX_ITER);
        
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
    iter = iter + 1;
    
    if SHOW_CONVERG
        fBestCurve(iter,1) = gBestFit;
        fesCurve(iter,1) = fes;
    end
end
xBest = gBest;
fBest = gBestFit;
end