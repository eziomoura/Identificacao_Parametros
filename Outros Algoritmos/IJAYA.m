function [xBest, fBest, fBestCurve, fesCurve] = IJAYA(fobj, LB, UB, PARAM, MAX_FES, SHOW_CONVERG)
% Descrição
%     IJAYA minimiza a fobj usando a metaheurística "Improved JAYA",
% conforme descrita em [1].
% Entradas:
%   fobj - Função objetivo a ser minimizada
%   LB - Vetor linha com os limites inferiores de cada parâmetro
%   UB - Vetor linha com os limites superior de cada parâmetro
%   PARAM - Estrutura com o seguintes campos:
%      pop - Tamanho da população
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
%   [1] YU, K.; LIANG, J. J.; QU, B. Y.; CHEN, X.; WANG, H. Parameters identification of photovoltaic models using an improved JAYA optimization algorithm. Energy Conversion and Management, v. 150, n. July, p. 742–753, 2017. 
%% Populacao inicial
POP_SIZE = PARAM.pop;
DIM = length(LB);
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
fit = fobj(x);     % Avalicao do fitness de cada individuo
fes = POP_SIZE;    % Quantidade de avaliações da função objetivo

%% Dados da curva de convergência
if SHOW_CONVERG
    MAX_ITER = floor((MAX_FES - fes)/(POP_SIZE)) + 1;
    fBestCurve = NaN(MAX_ITER,1);
    fesCurve = NaN(MAX_ITER,1);
    fBestCurve(1) = min(fit);
    fesCurve(1) = fes;
else
    fBestCurve = [];
    fesCurve = [];
end
z = rand; % Inicializacao do mapa logistico
iter = 1; % contador de iterações
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
    end % fim da atualiação a posicao de cada individuo
    iter = iter +1;
    if SHOW_CONVERG
        fBestCurve(iter,1) = min(fit);
        fesCurve(iter,1) = fes;
    end
end % encerra quantidade maxima de avaliacoes da funcao objetivo
% Remove NaNs
fBestCurve = fBestCurve(1:iter, 1);
fesCurve = fesCurve(1:iter, 1);
%
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
