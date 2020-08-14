function f = getAlgo(algoSelected)
switch algoSelected
    % DE
    case 'JADE'
        f = @JADE;
    case 'EJADE'
        f = @EJADE;
    case 'DE'
        f = @DE;
        % TLBO
    case 'TLBO'
        f = @TLBO;
    case 'ITLBO'
        f = @ITLBO;
        % BFS
    case 'BFS'
        f = @BFS;
    case 'ABC'
        f = @ABC;
    case 'PSO'
        f = @PSO;
    case 'IJAYA'
        f = @IJAYA;
    case 'PGJAYA'
        f = @PGJAYA;
    otherwise
        error('algoritmo não encontrado');
end