function [limite_inf, limite_sup] = getLimit(model, Ns, curva)
% ref: PGJAYA
% 1D = [Iph, I0, n, Rs, Rp]; 
name = curva.name;
switch model
    case '1D' % modelo um diodo
        if Ns > 1 % modulo
            switch name
                case {'Photowatt-PWP 201 '}
                    limite = [0 2 % Iph
                        0 50e-6   % I0
                        1 50/Ns   % n
                        0 2       % Rs
                        0 2000];  % Rp
                otherwise
                    limite = [0 max(curva.I)*1.05 %Iph
                        0 1e-6   %I01
                        1 2   %n1
                        0 0.3     %Rs
                        0 1000];  %Rp
            end  
            
        else % célula
            limite = [0 1
                0 1e-6
                1 2
                0 0.5
                0 100];
        end
   case '2D' % modelo dois diodos  {Iph, I01, I02, n1, n2, Rs, Rp}
       if Ns > 1 % modulo
           switch name
               case {'Photowatt-PWP 201 '}
                   limite = [0 2 %Iph
                       0 50e-6   %I01
                       0 50e-6   %I02
                       1 50/Ns   %n1
                       1 50/Ns   %n2
                       0 2       %Rs
                       0 2000];  %Rp
               otherwise
                   limite = [0 max(curva.I)*1.05 %Iph
                       0 1e-6   %I01
                       0 1e-6   %I02
                       1 50/Ns   %n1
                       1 50/Ns   %n2
                       0 0.3       %Rs
                       0 1000];  %Rp
           end
                    

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
