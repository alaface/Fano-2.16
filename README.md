# Automorphism Groups of Fano Threefolds of Family 2.16

This repository contains Magma code used to study the automorphism groups of Fano threefolds in family 2.16. These varieties are obtained by blowing up the complete intersection of two smooth quadrics in $\mathbb{P}^5$ along a conic. The scripts compute the possible group actions that preserve both the intersection and the conic, and identify the finite groups that can arise as automorphism groups of the resulting blow-ups.

## Contents

### `library.m`
A utility library with general-purpose functions used across all scripts:
- `FindLis(g, FF)`: Given a permutation $g$, returns a list of subgroups of $\mathrm{GL}_6(FF)$ preserving the matrix of absolute values and containing $-I$.
- `WedgePower`: Computes the $k$-th exterior power of a matrix.
- `Fiber`: Computes the fiber of a morphism over a given subscheme.

### `Grassmannian subvarieties.txt`
Sets up the ambient data for computations:
- The Grassmannian $\mathrm{Gr}(3,6) \subset \mathbb{P}^{19}$, via the Plücker embedding.
- A subvariety $W \subset \mathrm{Gr}(3,6)$ defined by Plücker equations and additional constraints.
- Coordinate rings, ambient maps, and projections for use in fixed point calculations.

### Case folders
Each folder corresponds to a specific conjugacy class of permutations in $S_6$ and contains a script that computes the fixed loci of planes under the associated subgroup action.

- `case(12)(34)(56)`  
- `case(123)(456)`  
- `case(12345)`  
- `case(123456)`  

Each script performs the following steps:
1. Initializes a cyclotomic or quadratic field suitable for the symmetry.
2. Uses `FindLis` to determine relevant subgroups of $\mathrm{GL}_6$.
3. Computes fixed loci of these subgroups acting on $W \subset \mathrm{Gr}(3,6)$.
4. Filters the fixed loci by non-emptiness.
5. Outputs the dimension of each fixed locus and identifies the group acting trivially on it.
6. For selected loci, constructs the corresponding fixed planes in $\mathbb{P}^5$.

### `Permutation scheme`
Magma script that constructs, for a given permutation $\rho$, the scheme of points $(a_0: \dots : a_5) \in \mathbb{P}^5$ such that the pencil of quadrics defined by two diagonal quadrics is invariant under $\rho$.

- The script builds a matrix involving the coefficients $a_i$ and their permutation under $\rho$, and imposes a rank 2 condition via minors.
- It excludes degenerate points where coefficients vanish or coincide.
- The result is a scheme parameterizing pencils of quadrics invariant under $\rho$.

## Usage

To run one of the case studies:

```magma
FF<i> := CyclotomicField(n); // or QuadraticField(-1)
load "library.m";
load "Grassmannian subvarieties.txt";

// Replace below with your permutation
g := Sym(6)!(1,2)(3,4)(5,6);
lis := FindLis(g, FF);

// Then run the corresponding script logic...
