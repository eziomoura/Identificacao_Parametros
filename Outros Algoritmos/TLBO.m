function [xBest, fBest, fBestCurve, fesCurve] = TLBO(fobj, LB, UB, PARAM, MAX_FES, seeConverg)
% Descrição
%     TLBO miniza a fobj usando a metaheurística teaching-learning-based optimization,
% conforme descrita em [1].
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
%   [1] VENKATA RAO, R. Review of applications of tlbo algorithm and a tutorial for beginners to solve the unconstrained and constrained optimization problems. Decision Science Letters, v. 5, n. 1, p. 1–30, 2016. 

%%
POP_SIZE = PARAM.pop;
DIM = length(LB);     % qtd de variaveis de design
% Inicializa a população
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
fit = fobj(x);

fes = POP_SIZE; % contador de avaliacoes da funcao objetivo

% pre-alocacao de memoria
MAX_ITER = floor((MAX_FES - fes)/(2*POP_SIZE)) + 1;

% pre alocacao da curva de convergência
if seeConverg
    fBestCurve = NaN(MAX_ITER, 1);
    fesCurve =  NaN(MAX_ITER, 1);
    fBestCurve(1) = min(fit);
    fesCurve(1) = fes;
end

iter = 1; % contador de iteracoes
while(fes + 2*POP_SIZE <= MAX_FES)
    for i = 1:POP_SIZE
       %% Teacher Phase
        xMean = mean(x);
        
        % determina o professor
        [~,id] = min(fit);
        xBest = x(id,:);
        
        % Teaching factor
        TF = randi([1 2]);
        
        % gera nova solucao
        xNew = x(i,:) + rand(1,DIM).*(xBest - TF*xMean);
        
        % verifica os limites
        xNew = boudaryCorrection(xNew, LB, UB, DIM, 1);        
        
        % avalia a nova solucao
        fitNew = fobj(xNew);
        if (fitNew < fit(i))
            x(i,:) = xNew;
            fit(i) = fitNew;
        end
      %% Learner Phase
        % seleciona um parceiro aleatoriamente
        randNums = randperm(POP_SIZE,2);
        if randNums(1) ~= i
            idRand = randNums(1);
        else
            idRand = randNums(2);
        end
        
        % gera nova solução
        if (fit(i) < fit(idRand))
            xNew = x(i,:) + rand(1, DIM).*(x(i,:)- x(idRand,:));
        else
            xNew = x(i,:) + rand(1, DIM).*(x(idRand,:)- x(i,:));
        end
        
        % verifica os limitantes
        xNew = boudaryCorrection(xNew, LB, UB, DIM, 1);
        
        % avalia a nova solução
        fitNew = fobj(xNew);
        if(fitNew <  fit(i))
            x(i,:) = xNew;
            fit(i) = fitNew;
        end        
    end
    
    fes = fes + 2*POP_SIZE;
    iter = iter + 1;
    
    if seeConverg
        fBestCurve(iter) = min(fit);
        fesCurve(iter) = fes;
    end  
end
% Remove os NaN em excesso
fBestCurve = fBestCurve(1:iter, 1);
fesCurve = fesCurve(1:iter, 1);

% retorna a melhor solução
[fBest, id] = min(fit);
xBest = x(id,:);
end
%%
function xNew = boudaryCorrection(xNew, lb, ub, DIM, popSize)
LBmatrix = repmat(lb, popSize,1);
UBmatrix = repmat(ub, popSize,1);

u = (xNew < LBmatrix);
xNew(u) = LBmatrix(u);
u = xNew > UBmatrix;
xNew(u) = UBmatrix(u);
end