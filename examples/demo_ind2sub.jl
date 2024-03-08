using Comodo

A = rand(30)
B = rand(5,6)
C = rand(3,5,2)

ind = [1,2,3,4,8,12,30]
IJK_A = ind2sub(size(A),ind)
IJK_B = ind2sub(size(B),ind)
IJK_C = ind2sub(size(C),ind)

l_A = all([A[ind[i]] == A[IJK_A[i][1]] for i ∈ eachindex(ind)])
l_B = all([B[ind[i]] == B[IJK_B[i][1],IJK_B[i][2]] for i ∈ eachindex(ind)])
l_C = all([C[ind[i]] == C[IJK_C[i][1],IJK_C[i][2],IJK_C[i][3]] for i ∈ eachindex(ind)])
