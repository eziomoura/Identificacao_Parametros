function [xBest, fBest, fBestCurve, fesCurve] = PSO(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
% Descri��o
%     XXXX miniza a fobj usando a metaheur�stica XXXXX,
% conforme descrita em [1] e [2].
% Entradas:
%   fobj - Fun��o objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada par�metro
%   UB - Vetor linha com os limites superior de cada par�metro
%   POP_SIZE - Inteiro com o tamanho da popula��o
%   MAX_FES - Inteiro com o quantidade m�xima de avali��es da fun��o objetivo
%   showConverg - Valor boleano que se for VERDADEIRO, ativar� as sa�das com os vetores 
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