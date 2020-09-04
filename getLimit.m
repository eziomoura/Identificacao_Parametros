function [limite_inf, limite_sup] = getLimit(funType, Ns)
%          1D = [Iph, I0, n, Rs, Rp]; 
% 2D {'Iph', 'I01', 'I02', 'n1','n2', 'Rs', 'Rp'}
switch funType
    case 1
        if Ns > 1
            limite_inf = [0, 0, 1, 0, 0];            % limite inferior do módulo
            limite_sup = [2, 50e-6, 50/Ns, 2, 2000]; % limite superior
        else
            limite_inf = [0, 0, 1, 0, 0];            % limite inferior da célula
            limite_sup = [1, 1e-6, 2, 0.5, 100];     % limite superior
        end
   case 2
        if Ns > 1
            limite_inf = [0, 0, 1, 0, 0];            % limite inferior do módulo
            limite_sup = [2, 50e-6, 50/Ns, 2, 2000]; % limite superior
        else
            limite_inf = [0, 0, 1, 0, 0];            % limite inferior da célula
            limite_sup = [1, 1e-6, 2, 0.5, 100];     % limite superior
        end
    otherwise
        error('limite não encontrado')
end
end