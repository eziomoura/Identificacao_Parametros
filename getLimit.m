function [limite_inf, limite_sup] = getLimit(model, Ns)
% ref: PGJAYA
% 1D = [Iph, I0, n, Rs, Rp]; 

switch model
    case '1D' % modelo um diodo
        if Ns > 1 % modulo
            limite = [0 2 % Iph
                0 50e-6   % I0
                1 50/Ns   % n
                0 2       % Rs
                0 2000];  % Rp

            
        else % célula
            limite = [0 1
                0 1e-6
                1 2
                0 0.5
                0 100];
        end
   case '2D' % modelo dois diodos  {Iph, I01, I02, n1, n2, Rs, Rp}
        if Ns > 1 % modulo
            limite = [0 2 %Iph
                0 50e-6   %I01
                0 50e-6   %I02
                1 50/Ns   %n1
                1 50/Ns   %n2
                0 2       %Rs
                0 2000];  %Rp

        else % celula 
            limite = [0 1
                0 1e-6
                0 1e-6
                1 2
                1 2
                0 0.5
                0 100];

        end
    otherwise
        error('limite não encontrado')
end
limite_inf(1,:) = limite(:,1);
limite_sup(1,:) = limite(:,2);
end
% STP6-120/36
% limite =[0 8
%     0 50
%     1 50
%     0 0.36
%     0 1500];

% STM6-40/36
% limite =[0 2
% 0 50e-6
% 1 60
% 0 0.36
% 0 1000];
