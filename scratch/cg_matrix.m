function CG = cg_matrix (l1 , l2)

% Initialize CG - matrix with zero entries .
n = (2* l1 +1) * (2* l2 +1);
CG = zeros (n,n);

% Fill CG - matrix band - wise ...
for l=abs(l2 -l1 ): l2+l1
    % ... and for an increasing value m ...
    for m=-l:l
        % ... with CG - coefficients .
        [rows , cols , cg] = cg_coefficients (l,m,l1 ,l2 );
        for i=1: length (cg)
            CG( rows (i), cols (i)) = cg(i);
        end
    end
end

end