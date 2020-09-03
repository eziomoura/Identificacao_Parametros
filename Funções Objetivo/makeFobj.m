classdef makeFobj
    properties
        Vmed
        Imed
        Ns
        T
        k = 1.3806503e-23;    % Boltzmann [J/K] 1.380649 [1.38065040000000e-23]
        q = 1.60217646e-19;   % Electron charge [C]
        Vt                    % tens�o t�rmica
        metrica
        modelo
        grandeza
    end
    
    properties (Access = private)
        f                     % fun��o objetivo
        typeFobj              % identificador da fun��o objetivo
    end
    
    methods
        function this = makeFobj(Vmed, Imed, Ns, T, typeFobj)
            this.Vmed = Vmed;
            this.Imed = Imed;
            this.Ns = Ns;
            this.T = T + 273.15;
            this.Vt = Ns * this.k * this.T / this.q;
            switch typeFobj
                case 1
                    this.f = @RMSE_CURRENT_SINGLE_DIODE;
                    this.metrica = 'RMSE';
                    this.modelo = '1D';
                    this.grandeza = 'CURRENT';
                case 2 
                    this.f = @RMSE_POWER_ONE_DIODE;
                    this.metrica = 'RMSE';
                    this.modelo = '1D';
                    this.grandeza = 'POWER';
                    
                otherwise
                    error('fun��o objetivo n�o encontrada');
            end
        end
        
        function y = Objective(this, x) 
            [POP_SIZE,~] = size(x);
             y = this.f(x, this, POP_SIZE);
        end
    end
end