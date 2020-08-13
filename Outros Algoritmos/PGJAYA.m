function [Iph, I0, n, Rs, Rp, RMSE, converg_RMSE, converg_fes] = PGJAYA(fobj, LB, UB, POP_SIZE, MAX_FES, seeConverg)
%% constantes
DIM = length(LB);
%% Populacao inicial
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
fit = fobj(x);     % Avalicao do fitness de cada individuo
fes = POP_SIZE; % Quantidade de avaliações da função objetivo

%% Dados da curva de convergência
if seeConverg
    MAX_ITER = floor((MAX_FES - POP_SIZE)/POP_SIZE);
    converg_RMSE = zeros(MAX_ITER +1,1);
    converg_fes = zeros(MAX_ITER +1,1);
    converg_RMSE(1) = sqrt(min(fit));
    converg_fes(1) = fes;
end
z = rand; % Inicializacao do mapa logistico
iter = 1; % contador de iterações

R = POP_SIZE - (1:POP_SIZE);
P = (R/POP_SIZE).^2;      % vetor de probabilidades
while(fes + POP_SIZE <= MAX_FES)
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
    end % fim da atualiação a posicao de cada individuo
    
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
    if seeConverg
        converg_RMSE(iter,1) = sqrt(min(fit));
        converg_fes(iter,1) = fes;
    end
end % encerra quantidade maxima de avaliacoes da funcao objetivo
[MSE, id] = min(fit);
RMSE = sqrt(MSE);
xBest = x(id,:);
Iph = xBest(1);
I0 = xBest(2);
n = xBest(3);
Rs = xBest(4);
Rp = xBest(5);
end

function xNew = checkBoundary(xNew, lb, ub, popSize, dim)
u = (xNew < lb.*ones(popSize,dim)) | (xNew > ub.*ones(popSize,dim));
randomMatrix = lb.*ones(popSize, dim) + (ub - lb).*ones(popSize, dim).*rand(popSize, dim);
xNew(u) = randomMatrix(u);
end
