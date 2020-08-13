%%
% o autor não indica em qual etapa faz a restrição das variáveis
%%
function [Iph, I0, n, Rs, Rp, RMSE, RMSE_curve, fes_curve] = ITLBO(fobj, LB, UB, POP_SIZE, MAX_FES, seeConverg)
% Vmed e Imed devem ser vetores coluna
% LB e UB devem ser vetores linha [1,dim]
DIM = length(LB);     % qtd de variaveis de design
%%
% Generation of initial population
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
fit = fobj(x);

fes = POP_SIZE; % contador de avaliacoes da funcao objetivo

% pre-alocacao de memoria
MAX_ITER = floor((MAX_FES - POP_SIZE)/(2*POP_SIZE));

% pre alocacao da curva de convergência
if seeConverg
    RMSE_curve = zeros(MAX_ITER + 1, 1);
    fes_curve =  zeros(MAX_ITER + 1, 1);
    RMSE_curve(1) = min(fit);
    fes_curve(1) = fes;
end

iter = 1; % contador de iteracoes
while(fes+2*POP_SIZE <= MAX_FES)
    for i = 1:POP_SIZE
        %% Teacher Phase
        xMean = mean(x);
        fitMean = fobj(xMean);
        
        % Determination of teacher
        [~,id] = min(fit);
        xBest = x(id,:);
        
        % Teaching factor
        TF = randi([1 2]);
        
        % Generation of a new solution
        if fit(i) < fitMean
            idRand = randperm(POP_SIZE, 3);
            idRand(idRand == i) = [];
            xNew = x(i,:) + rand*(xBest - x(i,:))...
                          + rand*(x(idRand(1),:) - x(idRand(2),:));
        else
            xNew = x(i,:) + rand*(xBest - TF*xMean);
        end
        
        % Bounding of the solution
        xNew = boudaryCorrection(xNew, LB, UB, DIM, 1);        
        
        % Evaluation of objective function
        fitNew = fobj(xNew);
        
        % Greedy selection
        if (fitNew < fit(i))
            x(i,:) = xNew;
            fit(i) = fitNew;
        end
%%      Learner Phase
        % selecionar parceiros aleatoriamente
        idRand = randperm(POP_SIZE,5);
        idRand(idRand==i) = [];
        
        % Generation of a new solution
        if (fit(i) < fitMean)
            xNew = x(i,:) + rand*(x(idRand(1),:) - x(idRand(2),:));
        else
            xNew = x(i,:) + rand*(x(idRand(1),:) - x(idRand(2),:))...
                          + rand*(x(idRand(3),:) - x(idRand(4),:));
        end
        
        % Bounding of the solution
        xNew = boudaryCorrection(xNew, LB, UB, DIM, 1);
        
        % Evaluation of objective function
        fitNew = fobj(xNew);
        
        % Greedy selection
        if(fitNew <  fit(i))
            x(i,:) = xNew;
            fit(i) = fitNew;
        end        
    end
    
    fes = fes + 2*POP_SIZE;
    iter = iter + 1;
    
    if seeConverg
        RMSE_curve(iter) = min(fit);
        fes_curve(iter) = fes;
    end  
end
% Extracting the best solution
[MSE, id] = min(fit);
RMSE = sqrt(MSE);
RMSE_curve = sqrt(RMSE_curve);
xBest = x(id,:);
Iph = xBest(1);
I0 = xBest(2);
n = xBest(3);
Rs = xBest(4);
Rp = xBest(5);
end
%%
function xNew = boudaryCorrection(xNew, LB, UB, DIM, POP_SIZE)
%% LB e UB devem ser matrizes com dimensao [POP_SIZE, dim]
u = (xNew < LB) | (xNew > UB);
randomMatrix = LB + (UB - LB).*rand(POP_SIZE, DIM);
xNew(u) = randomMatrix(u);
end
