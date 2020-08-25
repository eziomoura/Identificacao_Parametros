function [f, param] = getAlgo(algoSelected, paramData)
switch algoSelected
    case 'BFS'
        f = @BFS;
        param = paramData.BFS;

    case 'JADE'
        f = @JADE;
        param = paramData.JADE;
    case 'EJADE'
        f = @EJADE;
        param = paramData.EJADE;
    case 'DE'
        f = @DE;
        param = paramData.DE;

    case 'TLBO'
        f = @TLBO;
        param = paramData.TLBO;
    case 'ITLBO'
        f = @ITLBO;
        param = paramData.ITLBO;
    case 'TLABC'
        f = @TLABC;
        param = paramData.TLABC;

    case 'ABC'
        f = @ABC;
        param = paramData.ABC;
    case 'CIABC'
        f = @CIABC;
        param = paramData.CIABC;

    case 'PSO'
        f = @PSO;
        param = paramData.PSO;

    case 'IJAYA'
        f = @IJAYA;
        param = paramData.IJAYA;
    case 'PGJAYA'
        f = @PGJAYA;
        param = paramData.PGJAYA;
    otherwise
        error('algoritmo não encontrado');
end