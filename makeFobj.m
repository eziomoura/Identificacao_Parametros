classdef makeFobj
    properties
        Vmed
        Imed
        Ns
        T
        k = 1.3806503e-23;    % Boltzmann [J/K] 1.380649 [1.38065040000000e-23]
        q = 1.60217646e-19;   % Electron charge [C]
        Vt
    end
    
    methods
        function obj = makeFobj(Vmed, Imed, Ns, T)
            obj.Vmed = Vmed;
            obj.Imed = Imed;
            obj.Ns = Ns;
            obj.T = T + 273.15;
            obj.Vt = Ns*obj.k*obj.T/obj.q;
        end
        
        function y = Fobj(this, x) 
            [POP_SIZE,~] = size(x);
             y = RMSE_CURRENT_ONE_DIODE(x, this.Vmed, this.Imed, this.Vt, POP_SIZE);
        end
    end
end