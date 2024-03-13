using Comodo

A = rand(30)
B = rand(5,6)
C = rand(3,5,2)

ind = [1,2,3,4,8,12,30]

IJK_A = [[1],[2],[3],[4],[8],[12],[30]]
IJK_B = [[1, 1], [2, 1], [3, 1], [4, 1], [3, 2], [2, 3], [5, 6]]
IJK_C = [[1, 1, 1], [2, 1, 1], [3, 1, 1], [1, 2, 1], [2, 3, 1], [3, 4, 1], [3, 5, 2]]

all(sub2ind(size(A),IJK_A)==ind)
all(sub2ind(size(B),IJK_B)==ind)
all(sub2ind(size(C),IJK_C)==ind)