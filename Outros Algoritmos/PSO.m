% Developer: Ezio Moura
% PSO code based on "Parameter extraction of solar cells using particle swarm
% optimization" by Meiying Ye, Xiaodong Wang, and Yousheng Xu
%%
function [Iph, I0, n, Rs, Rp, RMSE, converg_RMSE, converg_fes] = PSO(fobj, LB, UB, POP_SIZE, MAX_FES, seeConverg)
%% parâmetros do algoritmo
dim = length(LB); % qtd de variaveis de design
c1 = 2;
c2 = 2;
vMax = (UB - LB);
w = 1;
%% inicialize positions and velocity
DIM = length(LB);
x = LB + (UB - LB).*rand(POP_SIZE, DIM);
v = zeros(size(x));
fobjValue = fobj(x);
pBest = x; % personal best of each particle
pBestFit = fobjValue;
[gBestFit, idBest] = min(fobjValue);
gBest = x(idBest,:); % global best

%% pre alocacao
fes = POP_SIZE;
if seeConverg
    MAX_ITER = floor((MAX_FES - POP_SIZE)/POP_SIZE);
    converg_RMSE = zeros(MAX_ITER + 1, 1);
    converg_fes = zeros(MAX_ITER + 1, 1);
    converg_RMSE(1) = sqrt(min(fobjValue));
    converg_fes(1) = fes;
end
iter = 1;
while(fes+POP_SIZE <= MAX_FES)
    for i = 1:POP_SIZE
        % update velocity
        v(i,:) = w*v(i,:) + c1*rand(1,dim).*(pBest(i,:) - x(i,:)) + c2*rand(1,dim).*(gBest - x(i,:));
        v(i,:) = min(v(i,:),  vMax);
        v(i,:) = max(v(i,:), -vMax);
        
        % update position
        x(i,:) = x(i,:) + v(i,:);
        x(i,:) = min(x(i,:), UB);
        x(i,:) = max(x(i,:), LB);
        
        % calc inercia weight
        w = 0.9 -(0.4)/(iter/MAX_ITER);
        
        % update personal best
        fobjValue(i) = fobj(x(i,:));
        if fobjValue(i) < pBestFit(i)
            pBest(i,:) = x(i,:);
            pBestFit(i) = fobjValue(i);
            
            % update global best (so far)
            if fobjValue(i) < gBestFit
                gBest = x(i,:);
                gBestFit = fobjValue(i);
            end
        end
    end
    
    fes = fes + POP_SIZE;
    iter = iter + 1;
    
    if seeConverg
        converg_RMSE(iter,1) = sqrt(gBestFit);
        converg_fes(iter,1) = fes;
    end
end

RMSE = sqrt(gBestFit);
Iph = gBest(1);
I0 = gBest(2);
n = gBest(3);
Rs = gBest(4);
Rp = gBest(5);
end