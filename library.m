// Computes the intersection product of two divisors using an intersection matrix
// Input: vectors a, b in a vector space; symmetric matrix M
// Output: scalar value of the intersection product
qua := function(a,b,M)
    K := BaseRing(Parent(M));
    u := Eltseq(a);
    v := Eltseq(b);
    return (Matrix(K,1,#u,u)*M*Matrix(K,#v,1,v))[1,1];
end function;

// Returns the matrix of absolute values of the entries of M
// Input: matrix M over a finite field FF
// Output: matrix with entries |M_{i,j}| in FF
AbsMatrix := function(M,FF)
    return Matrix([[Norm(FF!M[i,j]) : j in [1..Ncols(M)]] : i in [1..Nrows(M)]]);
end function;

// Finds subgroups of H whose elements have absolute value matrices equal to M
// and that contain -I (minus the identity matrix)
// Input: permutation g generating a subgroup of GL(n, FF)
// Output: list of subgroups satisfying the condition
FindLis := function(g)
    G := Parent(g);
    M := PermutationMatrix(FF, g);
    Glin := GL(Nrows(M), FF);

    // Generate diagonal sign-change matrices
    Dgens := [
        DiagonalMatrix(FF, [(j eq i select -1 else 1) : j in [1..6]])
        : i in [1..6]
    ];

    H := sub< Glin | Dgens, M >;

    // Extract subgroups of H
    lis := [U`subgroup : U in Subgroups(H)];

    // Keep only those containing elements with absolute value matrix equal to M
    lis := [U : U in lis | &or[AbsMatrix(s) eq M : s in U]];

    I := IdentityMatrix(FF, 6);

    // Return only those subgroups containing -I
    return [U : U in lis | -I in U];
end function;

// Computes the k-th exterior power of a square matrix P
// Input: matrix P and integer k
// Output: matrix representing the k-th exterior power
WedgePower := function(P, k)
    K := BaseRing(P);
    n := NumberOfRows(P);
    indices := Subsets({1..n}, k);
    Wedge_P := ZeroMatrix(K, #indices, #indices);
    index_list := Sort([Setseq(u) : u in indices]);

    for i in [1..#index_list] do
        for j in [1..#index_list] do
            rows := index_list[i];
            cols := index_list[j];
            Wedge_P[i, j] := Determinant(Submatrix(P, rows, cols));
        end for;
    end for;

    return Wedge_P;
end function;

// Computes the fiber of a morphism f over a scheme Z
// Input: morphism f: X -> Y, scheme Z in Y
// Output: scheme in X representing the fiber over Z
Fiber := function(f,Z)
    equ := DefiningEquations(f);
    bas := MinimalBasis(Z);
    P := Domain(f);
    X := Scheme(P,[Evaluate(g,equ) : g in bas]);
    return Complement(X,Scheme(P,equ));
end function;
