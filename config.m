Param.ABC.pop = 50;
Param.ABC.limit = 200;

Param.CIABC.pop = 50;
Param.CIABC.limit = 200;

Param.DE.pop = 50;
Param.DE.F = 0.95; % mutation factor
Param.DE.CR = 0.8; % crossover probability

Param.JADE.pop = 50; % ref: EJADE.m
Param.JADE.p = 0.05; % ref JADE - top p% melhores, determines the greediness of the mutation strategy
Param.JADE.c = 1/20; % ref JADE - controls the rate of parameter adaptation

Param.EJADE.popMax = 50; % ref: EJADE.m
Param.EJADE.popMin = 4;  % ref EJADE.m
Param.EJADE.p = 0.05; % ref JADE - top p% melhores, determines the greediness of the mutation strategy
Param.EJADE.c = 0.1; % ref JADE - controls the rate of parameter adaptation

Param.IJAYA.pop = 50;
Param.PGJAYA.pop = 50;

Param.ITLBO.pop = 50;
Param.TLBO.pop = 50;

Param.TLABC.pop = 50;
Param.TLABC.limit = 200;
Param.TLABC.F = rand;

Param.PSO.pop = 50;
Param.PSO.c1 = 2;
Param.PSO.c2 = 2;
Param.PSO.w = 1;
%Param.PSO.vmax = 50;


Param.BFS.pop = 50;