%%
% o autor n�o indica em qual etapa faz a restri��o das vari�veis
%%
function [xBest, fBest, fBestCurve, fesCurve] = ITLBO(fobj, LB, UB, POP_SIZE, MAX_FES, seeConverg)
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
DIM = length(LB);     % qtd de variaveis de design
%%
% Generation of initial population
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
        
        % Checar limites
        xNew = boudaryCorrection(xNew, LB, UB, DIM, 1);        
        
        % avaliar a nova solucao
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
        
        % Checar limites
        xNew = boudaryCorrection(xNew, LB, UB, DIM, 1);
        
        % Avaliar fun�a� objetivo
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
        fBestCurve(iter) = min(fit);
        fesCurve(iter) = fes;
    end  
end
% Melhor solu��o
[fBest, id] = min(fit);
xBest = x(id,:);
end
%%
function xNew = boudaryCorrection(xNew, LB, UB, DIM, POP_SIZE)
%% LB e UB devem ser matrizes com dimensao [POP_SIZE, dim]
u = (xNew < LB) | (xNew > UB);
randomMatrix = LB + (UB - LB).*rand(POP_SIZE, DIM);
xNew(u) = randomMatrix(u);
end
