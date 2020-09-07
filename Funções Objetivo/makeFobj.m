classdef makeFobj
    properties
        Vmed
        Imed
        Ns
        T
        k = 1.3806503e-23;    % Boltzmann [J/K] 1.380649 [1.38065040000000e-23]
        q = 1.60217646e-19;   % Electron charge [C]
        Vt                    % tensão térmica
        metrica
        modelo
        grandeza
    end
    
    properties (Access = private)
        f                     % função objetivo
        typeFobj              % identificador da função objetivo
    end
    
    methods
        function this = makeFobj(Vmed, Imed, Ns, T, obj_code)
            this.Vmed = Vmed;
            this.Imed = Imed;
            this.Ns = Ns;
            this.T = T + 273.15;
            this.Vt = Ns * this.k * this.T / this.q;
            this.modelo   = obj_code.modelo;
            this.metrica  = obj_code.metrica;
            this.grandeza = obj_code.grandeza;
            switch obj_code.modelo
                case '1D'
                    if (obj_code.metrica == 'RMSE' & obj_code.grandeza == 'I')
                        this.f = @RMSE_CURRENT_SINGLE_DIODE;
                    elseif(obj_code.metrica == 'RMSE' & obj_code.grandeza == 'P')
                        this.f = @RMSE_POWER_ONE_DIODE;
                    else
                        error('função objetivo não encontrada');
                    end
                        
                case '2D'
                    if (obj_code.metrica == 'RMSE' & obj_code.grandeza == 'I')
                        this.f = @RMSE_CURRENT_DOUBLE_DIODE;
                    end
                otherwise
                    error('função objetivo não encontrada');
            end
        end
        
        function y = Objective(this, x) 
            [POP_SIZE,~] = size(x);
             y = this.f(x, this, POP_SIZE);
        end
    end
end