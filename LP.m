function  Pt = LP(Adj, r, P0) 
deg_rowvec = ( sum(Adj,1) ).^0.5;  
deg_colvec = ( sum(Adj,2) ).^0.5;   
WAdj = (Adj./(deg_colvec+eps))./(deg_rowvec+eps) ;      
    Pt = P0;
    for T = 1:100 
        Pt1 = (1-r)*WAdj*Pt + r*P0;
        if norm(Pt1-Pt,1) < 10^-6
            break;
        end
        Pt = Pt1;
    end
end