function [xBest, fBest, fBest_curve, fes_curve] = IJAYA(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
% Descri��o
%     IJAYA minimiza a fobj usando a metaheur�stica "Improved JAYA",
% conforme descrita em [1].
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
%   [1] YU, K.; LIANG, J. J.; QU, B. Y.; CHEN, X.; WANG, H. Parameters identification of photovoltaic models using an improved JAYA optimization algorithm. Energy Conversion and Management, v. 150, n. July, p. 742�753, 2017. 
%% Populacao inicial
POP_SIZE = PARAM.pop
DIM = length(LB);
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
fit = fobj(x);     % Avalicao do fitness de cada individuo
fes = POP_SIZE;    % Quantidade de avalia��es da fun��o objetivo

%% Dados da curva de converg�ncia
if SHOW_CONVERG
    MAX_ITER = floor((MAX_FES - POP_SIZE)/POP_SIZE);
    fBest_curve = zeros(MAX_ITER +1,1);
    fes_curve = zeros(MAX_ITER +1,1);
    fBest_curve(1) = min(fit);
    fes_curve(1) = fes;
end
z = rand; % Inicializacao do mapa logistico
iter = 1; % contador de itera��es
while(fes + POP_SIZE <= MAX_FES)
    %% identifica o melhor e o pior
    [~,id] = sort(fit);
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
        % se nao for o melhor
        if i ~= id(1) 
            if rand < rand % usa eq 12: Self-adaptive weight
                xNew = x(i,:) + rand(1,DIM).*(xBest - abs(x(i,:))) - w*rand(1,DIM).*(xWorst - abs(x(i,:)));
            else % usa eq 13: Experience-based learning strategy
                randInd = randperm(POP_SIZE,3);
                randInd(i == randInd) = [];
                if fit(randInd(1)) < fit(randInd(2))
                    xNew = x(i,:) + rand(1,DIM).*(x(randInd(1),:) - x(randInd(2),:));
                else
                    xNew = x(i,:) - rand(1,DIM).*(x(randInd(1),:) - x(randInd(2),:));
                end
            end
        % se for o melhor
        else 
            z = 4*z*(1 - z);
            xNew = xBest + (2*z - 1).*rand(1,DIM);
        end
        
        % Faz checagem de bordas
        xNew = checkBoundary(xNew, LB, UB, 1, DIM);
        
        % avalia o novo individuo
        fitNew = fobj(xNew);
        fes = fes+1;
        if fitNew < fit(i)
            fit(i) = fitNew;
            x(i,:) = xNew;
        end
    end % fim da atualia��o a posicao de cada individuo
    iter = iter +1;
    if SHOW_CONVERG
        fBest_curve(iter,1) = min(fit);
        fes_curve(iter,1) = fes;
    end
end % encerra quantidade maxima de avaliacoes da funcao objetivo
[fBest, id] = min(fit);
xBest = x(id,:);
end

function xNew = checkBoundary(xNew, lb, ub, popSize, dim)
u = (xNew < lb.*ones(popSize,dim)) | (xNew > ub.*ones(popSize,dim));
randomMatrix = lb.*ones(popSize, dim) + (ub - lb).*ones(popSize, dim).*rand(popSize, dim);
xNew(u) = randomMatrix(u);
end
