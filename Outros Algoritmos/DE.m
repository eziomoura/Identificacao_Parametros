%DE
% Developer: Ezio Moura
%%
function [Iph, I0, n, Rs, Rp, RMSE, converg_RMSE, converg_fes] = DE(fobj, LB, UB, POP_SIZE, MAX_FES, seeConverg)
% Vmed e Imed devem ser vetores coluna
% LB e UB devem ser vetores linha [1,dim]
%% parâmetros do algoritmo
DIM = length(LB); % qtd de variaveis de design
F = 0.95;
CR = 0.8;
%maxIter = max_fes/(N); % numero maximo de iterações
%% Inicializa a populacao
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
fit = fobj(x);
fes = POP_SIZE;
if seeConverg
    MAX_ITER = (MAX_FES - POP_SIZE)/POP_SIZE;
    converg_RMSE = zeros(MAX_ITER,1);
    converg_fes = zeros(MAX_ITER,1);
    converg_RMSE(1) = sqrt(min(fit));
    converg_fes(1) = fes;
end
iter = 1;
while(fes + POP_SIZE <= MAX_FES)
    %% 1 - Mutation (DE/RAND/1)
    for i = 1:POP_SIZE
        r = randperm(POP_SIZE, 3);
        r = r(r ~= i);
        v(i,:) = x(i,:) + F.*(x(r(1),:) - x(r(2),:));
    end
    %% 2 - Crossover
    for i = 1:POP_SIZE
        p = randi(DIM);
        for j = 1:DIM
            if rand <= CR | j == p
                u(i,j) = v(i,j);
            else
                u(i,j) = x(i,j);
            end
        end
    end
    %% 
    % verificar limites
    xNew = u;
    xNew = boudaryCorrection(xNew, LB, UB, DIM, POP_SIZE);
    
    % avalia a nova posicao
    fitNew = fobj(xNew); % avalia nova posicao
    [x, fit] =  updatePosition(x, fit, xNew, fitNew);
    
    fes = fes + POP_SIZE;
    iter = iter +1;
    
    if seeConverg
        converg_RMSE(iter) = sqrt(min(fit));
        converg_fes(iter) = fes;
    end
end
[MSE, id] = min(fit);
RMSE = sqrt(MSE);
xBest = x(id,:);
Iph = xBest(1);
I0 = xBest(2);
n = xBest(3);
Rs = xBest(4);
Rp = xBest(5);
end

function xNew = boudaryCorrection(xNew, LB, UB, DIM, POP_SIZE)
%% LB e UB devem ser matrizes com dimensao [nBirds, dim]
    u = (xNew < LB) | (xNew > UB);
    randomMatrix = LB + (UB - LB).*rand(POP_SIZE, DIM);
    xNew(u) = randomMatrix(u);
end

function [x, fit] =  updatePosition(x, fit, xNew, fitNew)
    isBetter = fitNew < fit;
    if any(isBetter)
        fit(isBetter) = fitNew(isBetter);
        x(isBetter, :) = xNew(isBetter, :);
    end
end
