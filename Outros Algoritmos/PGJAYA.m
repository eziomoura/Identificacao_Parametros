function [xBest, fBest, converg_RMSE, converg_fes] = PGJAYA(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
% Descri��o
%     PGJAYA miniza a fobj usando a metaheur�stica "performance-guided JAYA",
% conforme descrita em [1]. 
%
% Entradas:
%   fobj - Fun��o objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada par�metro
%   UB - Vetor linha com os limites superior de cada par�metro
%   PARAM - Estrutura com o seguintes campos:
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
%   [1] YU, K.; QU, B.; YUE, C.; GE, S.; CHEN, X.; LIANG, J. A performance-guided JAYA algorithm for parameters identification of photovoltaic cell and module. Applied Energy, v. 237, n. December 2018, p. 241�257, 2019. 

%% Populacao inicial
POP_SIZE = PARAM.pop;
DIM = length(LB); % dimens�o do problema
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
fit = fobj(x);     % Avalicao do fitness de cada individuo
fes = POP_SIZE; % Quantidade de avalia��es da fun��o objetivo

%% Dados para a curva de converg�ncia
if SHOW_CONVERG
    MAX_ITER = floor((MAX_FES - fes)/POP_SIZE) + 1;
    converg_RMSE = NaN(MAX_ITER, 1);
    converg_fes = NaN(MAX_ITER, 1);
    converg_RMSE(1) = min(fit);
    converg_fes(1) = fes;
else
    converg_RMSE(1) = NaN;
    converg_fes(1) = NaN;   
end
z = rand; % Inicializacao do mapa logistico
iter = 1; % contador de itera��es

R = POP_SIZE - (1:POP_SIZE);
P = (R/POP_SIZE).^2;      % vetor de probabilidades
while(fes + POP_SIZE +1 <= MAX_FES)
    %% identifica o melhor e o pior
    [~,id] = sort(fit);
    prob(id) = P;
    xBest = x(id(1),:);
    xWorst = x(id(end),:);
    %% funcao peso
    if(fit(id(end)) == 0)
        w = 1;
    else
        w = abs(fit(id(1))/(fit(id(end))))^2;
    end
    %% atualiza a posicao de cada individuo
    for i = 1:POP_SIZE
        if prob(i) < rand
            xNew = x(i,:) + rand(1,DIM).*(xBest - abs(x(i,:))) - w*rand(1,DIM).*(xWorst - abs(x(i,:)));
        else
            randInd(1) = randi(POP_SIZE);
            while randInd(1) == i || rand > prob(randInd(1))
                randInd(1)= randi(POP_SIZE);
            end
            randInd(2) = randi(POP_SIZE);
            while randInd(2)==i || randInd(2)==randInd(1)
                randInd(2)= randi(POP_SIZE);
            end
            xNew = x(i,:) + rand(1,DIM).*(x(randInd(1),:) - x(randInd(2),:));
        end
        % verifica limites
        xNew = checkBoundary(xNew, LB, UB, 1, DIM);
        
        % avalia a nova populacao
        fitNew = fobj(xNew);
        fes = fes+1;
        if fitNew < fit(i)
            fit(i) = fitNew;
            x(i,:) = xNew;
        end
    end % fim da atualia��o a posicao de cada individuo
    
    [~, id] = sort(fit);
    xBest = x(id(1),:);
    z = 4*z*(1 - z);
    for k = 1:DIM
        if rand < 1 - fes/MAX_FES
            xNew(k) = xBest(k) + rand*(2*z-1);
        else
            xNew(k) = xBest(k);
        end
    end
    
    % verifica limites
    xNew = checkBoundary(xNew, LB, UB, 1, DIM);
    
    % avalia a nova populacao
    fitNew = fobj(xNew);
    fes = fes+1;
    if fitNew < fit(id(end))
        fit(id(end)) = fitNew;
        x(id(end),:) = xNew;
    end
    
    iter = iter +1;
    if SHOW_CONVERG
        converg_RMSE(iter, 1) = min(fit);
        converg_fes(iter, 1) = fes;
    end
end % encerra quantidade maxima de avaliacoes da funcao objetivo
% Remove os NaN em excesso
converg_RMSE = converg_RMSE(1:iter, 1);
converg_fes = converg_fes(1:iter, 1);

[fBest, id] = min(fit);
xBest = x(id,:);
end

function xNew = checkBoundary(xNew, lb, ub, popSize, dim)
LBmatrix = repmat(lb, popSize,1);
UBmatrix = repmat(ub, popSize,1);

u = (xNew < LBmatrix);
xNew(u) = LBmatrix(u);
u = xNew > UBmatrix;
xNew(u) = UBmatrix(u);
end
