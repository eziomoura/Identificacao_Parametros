function [f, param] = getAlgo(algoSelected, paramData)
switch algoSelected
    case 'BFS'
        f = @BFS;
        param = paramData.BFS;
    case 'BFSmod'
        f = @BFSmod;
        param = paramData.BFSmod;

    case 'BFSmod3'
        f = @BFSmod3;
        param = paramData.BFSmod;

    case 'BFS3'
        f = @BFS3;
        param = paramData.BFSnew;
    case 'BFS4'
        f = @BFS4;
        param = paramData.BFSnew;
    case 'BFS5'
        f = @BFS5;
        param = paramData.BFSnew;
    case 'BFSnew'
        f = @BFSnew;
        param = paramData.BFSnew;
        
    case 'BFSnewBeta'
        f = @BFSnewBeta;
        param = paramData.BFSnew;
        
            case 'BFS2'
        f = @BFS2;
        param = paramData.BFSnew;
  %-------------      
    case 'SHADE'
        f = @SHADE;
        param = paramData.SHADE;  
    case 'MADE'
        f = @MADE;
        param = paramData.MADE;
    case 'SEDE'
        f = @SEDE;
        param = paramData.SEDE;
    case 'EJADE'
        f = @EJADE;
        param = paramData.EJADE;
    case 'JADE'
        f = @JADE;
        param = paramData.JADE;
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
    case 'ELPSO'
        f = @ELPSO;
        param = paramData.ELPSO;
        
    case 'IJAYA'
        f = @IJAYA;
        param = paramData.IJAYA;
    case 'PGJAYA'
        f = @PGJAYA;
        param = paramData.PGJAYA;
    otherwise
        error('algoritmo não encontrado');
end