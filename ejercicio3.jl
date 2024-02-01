#-----Patrick Townsend-----#
include("LUdescomp.jl")
A = [1 1 1;
     2 4 8;
     1 4 12]
b = [1,2,0]

L,U = LUfactor(A)
x = LU_linear_solve(L, U, b)
x2 = A\b

display(x)
display(x2)