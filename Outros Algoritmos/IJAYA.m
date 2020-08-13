function [Iph, I0, n, Rs, Rp, RMSE, converg_RMSE, converg_fes] = IJAYA(fobj, LB, UB, POP_SIZE, MAX_FES, seeConverg)
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
while(fes + POP_SIZE <= MAX_FES)
    %% identifica o melhor e o pior
    [~,id] = sort(fit);
    xBest = x(id(1),:);
    xWorst = x(id(end),:);
    %% funcao peso
    if(fit(id(end)) == 0)
        ww = 1;
    else
        ww = abs(fit(id(1))/(fit(id(end))))^2;
    end
    %% atualiza a posicao de cada individuo
    for i = 1:POP_SIZE 
        if i ~= id(1) % se nao for o melhor
            if rand < rand % usa eq 12: Self-adaptive weight
                xNew = x(i,:) + rand(1,DIM).*(xBest - abs(x(i,:))) - ww*rand(1,DIM).*(xWorst - abs(x(i,:)));
            else % usa eq 13: Experience-based learning strategy
                randInd = randperm(POP_SIZE,3);
                randInd(i == randInd) = [];
                if fit(randInd(1)) < fit(randInd(2))
                    xNew = x(i,:) + rand(1,DIM).*(x(randInd(1),:) -x(randInd(2),:));
                else
                    xNew = x(i,:) - rand(1,DIM).*(x(randInd(1),:) -x(randInd(2),:));
                end
            end
        else % se for o melhor
            z = 4*z*(1 - z);
            xNew = xBest + (2*z - 1).*rand(1,DIM);
        end
        
        % faz checagem de bordas
        xNew = checkBoundary(xNew, LB, UB, 1, DIM);
        
        % avalia a nova populacao
        fitNew = fobj(xNew);
        fes = fes+1;
        if fitNew < fit(i)
            fit(i) = fitNew;
            x(i,:) = xNew;
        end
    end % fim da atualiação a posicao de cada individuo
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
