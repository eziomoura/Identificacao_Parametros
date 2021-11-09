Param.ABC.pop = 50;
Param.ABC.limit = 200;

Param.CIABC.pop = 50;
Param.CIABC.limit = 200;

Param.DE.pop = 50;
Param.DE.F = 0.95; % mutation factor
Param.DE.CR = 0.8; % crossover probability

Param.SEDE.pop = 50;

Param.JADE.pop = 50; % ref: EJADE.m
Param.JADE.p = 0.05; % ref JADE - top p% melhores, determines the greediness of the mutation strategy
Param.JADE.c = 1/20; % ref JADE - controls the rate of parameter adaptation

Param.EJADE.popMax = 50; % ref: EJADE.m
Param.EJADE.popMin = 4;  % ref EJADE.m
Param.EJADE.p = 0.05; % ref JADE - top p% melhores, determines the greediness of the mutation strategy
Param.EJADE.c = 0.1; % ref JADE - controls the rate of parameter adaptation

Param.SHADE.pop = 50; % ref[1] em SHADE.m
Param.SHADE.H = 100;  % ref[1] em SHADE.m

Param.MADE.pop = 50; % ref[1] em MADE.m
Param.MADE.H = 100;  % ref[1] em MADE.m
Param.MADE.epsilon = 0.05;  % ref[1] em MADE.m, criterio para aplicar NM: f(x) < epsilon
%Param.MADE.epsilon = 0.01;  % para DDM



Param.IJAYA.pop = 50;
Param.PGJAYA.pop = 50;

Param.ITLBO.pop = 50;
Param.TLBO.pop = 50;

Param.TLABC.pop = 50;
Param.TLABC.limit = 200;
Param.TLABC.F = rand;

Param.PSO.pop = 60;
Param.PSO.c1 = 2;
Param.PSO.c2 = 2;
Param.PSO.w = 0.9;
%Param.PSO.vmax = 50;

Param.ELPSO.pop = 50; % ref[2] - ELPSO.m
Param.ELPSO.c1 = 2;  % personal acceleration coefficient
Param.ELPSO.c2 = 2;  % social acceleration coefficien
Param.ELPSO.w = 0.9; % inertia weight (decresse from 0.9 to 0.4)
Param.ELPSO.h = 1;   % standard deviation of Gaussian mutation
Param.ELPSO.s = 2;   % scale factor of cauchy distribuition
Param.ELPSO.F = 1.2; % scale factor of DE-based mutation


Param.BFS.pop = 50;
Param.BFSmod.pop = 50;
Param.BFSnew.pop = 50;
Param.BFSnew.popMax = 50; % 
Param.BFSnew.popMin = 4;  % 