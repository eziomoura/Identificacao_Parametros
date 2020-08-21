function f = getAlgo(algoSelected)
switch algoSelected
    case 'BFS'
        f = @BFS;
case 'JADEB'
        f = @JADEB;
    case 'JADE'
        f = @JADE;
    case 'EJADE'
        f = @EJADE;
    case 'DE'
        f = @DE;

    case 'TLBO'
        f = @TLBO;
    case 'ITLBO'
        f = @ITLBO;
    case 'TLABC'
        f = @TLABC;

    case 'ABC'
        f = @ABC;
    case 'CIABC'
        f = @CIABC;

    case 'PSO'
        f = @PSO;

    case 'IJAYA'
        f = @IJAYA;
    case 'PGJAYA'
        f = @PGJAYA;
    otherwise
        error('algoritmo não encontrado');
end