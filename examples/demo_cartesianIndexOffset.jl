using Comodo

#=
This demo shows the use of the `cartesianIndexOffset` function to create offsets
for CartesianIndex values. Essentially a CartesianIndex is created that can be 
added to an existing CartesianIndex to shift it. 
=#

# Offset in single direction
siz = [3, 4, 5] # Size of array 
d = 2 # Dimension for offset
offset = 2 # Offset amount
cn = cartesianIndexOffset(siz, d, offset) # The Cartesian index offset
println(cn)

# Offset in multiple directions
siz = [30, 40, 50] # Size of array 
d = [1,2,3] # Dimension for offset
offset = [2,-1,5] # Offset amount
cn = cartesianIndexOffset(siz, d, offset) # The Cartesian index offset
println(cn)