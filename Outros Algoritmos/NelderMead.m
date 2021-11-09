function [xbest,fbest, fBestCurve, fesCurve] = NelderMead(fobj, x0, tol, maxFes, SHOW_CONVERG)
% Descrição
%     NelderMead minimiza a fobj usando o método "Nelder-Mead Simplex Procedure"
% conforme descrito em [1]. Inicialização inspirada no metodo inidicado em
% [2].
%
% Entradas:
%   fobj - Função objetivo a ser minimizada
%   x0 - estimativa inicial
%   tol - tolerancia no valor minimo de fobj
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
%   [1] LAGARIAS, JEFFREY, C.; REEDS, JAMES, A.; WRIGHT, MARGARET, H.; WRIGHT, PAUL, E. Convergence properties of the nelder–mead simplex method in low dimensions. SIAM Journal on Optimization, v. 9, n. 1, p. 112–147, 1998. 
%   [2] http://var.scholarpedia.org/article/Nelder-Mead_algorithm
n = length(x0);
% Parametros
rho = 1;
chi = 2;
gama = 0.5;
sigma = 0.5;


% Inicializa os n+1 vértices com base em x0
s  = zeros(n+1,n); % s contem os n+1 vértices
s(1,:) = x0;

f = zeros(n+1,1);
f(1) = fobj(x0);

% inicializao dos outros vertices do simplex
step = 0.1*rand;
for i = 1:n
    s(i+1,:) = x0;
    if x0(i) == 0 %se for vetor nulo, não vai sair do lugar
        s(i+1,i) = step*10^-10; % coloca um valor muito pequeno
    else
        s(i+1,i) = x0(i) + step*x0(i);
    end
    f(i+1) = fobj(s(i+1,:));
end
fes = n+1;
iter = 1;
shrink = false;

% 1- Order
[f, id] = sort(f);
s = s(id,:);

if SHOW_CONVERG
    fBestCurve(iter,1) = f(1);
    fesCurve(iter,1) = fes;
end
while fes + (2+n) <= maxFes
    % 2- Reflect
    xcen = mean( s(1:n,:) ); % centroide dos n melhores
    xr = (1 + rho)*xcen - rho*s(end,:);
    fr = fobj(xr);
    fes = fes+1;
    if (f(1)<= fr & fr < f(n))
        s(end,:) = xr;
        f(end) = fr;
    else
        if fr < f(1)
            % 3- Expand
            xe = (1 + rho*chi)*xcen - rho*chi*s(end,:);
            fxe = fobj(xe);
            fes = fes + 1;
            if fxe < fr
                s(end,:) = xe;
                f(end) = fxe;
            else
                s(end,:) = xr;
                f(end) = fr;
            end
        else %(fr >= f(n)
            % 4 - contract
            if  fr < f(end)
                % a. Outside
                xc = (1 + rho*gama)*xcen - rho*gama*s(end,:);
                fc = fobj(xc);
                fes = fes+1;
                
                if fc <= fr
                    s(end,:) = xc;
                    f(end) = fc;
                else
                    shrink = true;
                end
            else
                % b. inside
                xcc = (1-gama)*xcen + gama*s(end,:);
                fcc = fobj(xcc);
                fes = fes+1;
                
                if fcc < f(end)
                    s(end,:) = xcc;
                    f(end) = fcc;
                else
                    shrink = true;
                end
            end
            if shrink
                for i = 2:n+1
                    s(i,:) = s(1,:) + sigma*(s(i,:) - s(1,:));
                    f(i) = fobj(s(i,:));
                end
                fes = fes + n;
            end
        end    
    end
    
    % 1- Order (para a proxima iteração)
    [f, id] = sort(f);
    s = s(id,:);
    
    % curva de convergencia
    iter = iter +1;
    if SHOW_CONVERG
        fBestCurve(iter,1) = f(1);
        fesCurve(iter,1) = fes;
    end
    
    % critério de parada
    if f(1) < tol
        break
    end
    shrink = false;
end
fbest = f(1);
xbest = s(1,:);
end