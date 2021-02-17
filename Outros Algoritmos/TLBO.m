function [xBest, fBest, fBestCurve, fesCurve] = TLBO(fobj, LB, UB, PARAM, MAX_FES, seeConverg)
% Descri��o
%     TLBO miniza a fobj usando a metaheur�stica XXXXX,
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

%%
POP_SIZE = PARAM.pop;
DIM = length(LB);     % qtd de variaveis de design
% Inicializa a popula��o
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
fit = fobj(x);

fes = POP_SIZE; % contador de avaliacoes da funcao objetivo

% pre-alocacao de memoria
MAX_ITER = floor((MAX_FES - POP_SIZE)/(2*POP_SIZE));

% pre alocacao da curva de converg�ncia
if seeConverg
    fBestCurve = zeros(MAX_ITER + 1, 1);
    fesCurve =  zeros(MAX_ITER + 1, 1);
    fBestCurve(1) = min(fit);
    fesCurve(1) = fes;
end

iter = 1; % contador de iteracoes
while(fes+2*POP_SIZE <= MAX_FES)
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
        
        % gera nova solu��o
        if (fit(i) < fit(idRand))
            xNew = x(i,:) + rand(1, DIM).*(x(i,:)- x(idRand,:));
        else
            xNew = x(i,:) + rand(1, DIM).*(x(idRand,:)- x(i,:));
        end
        
        % verifica os limitantes
        xNew = boudaryCorrection(xNew, LB, UB, DIM, 1);
        
        % avalia a nova solu��o
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
% retorna a melhor solu��o
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