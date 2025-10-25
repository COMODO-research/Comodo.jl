using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra

GLMakie.closeall()

# Visualization
fig = Figure(size=(1400,1000))
cartInd = CartesianIndices((3,4))
s = [0.5, 0.5, 0.5, 0.8, 0.5, 1.5, 1.0, 1.0, 0.5]
for (i_plot, tpms_type) in enumerate( (:P, :D, :D2, :N, :G, :IWP, :FRD, :S, :HG) )       
    F, V = tpms_sheet(tpms_type, s[i_plot]; x=range(0, 2*pi, 50), cap=true)        

    ij = cartInd[i_plot]
    ax1 = AxisGeom(fig[ij[1], ij[2]]; title=string(tpms_type))
    hp1 = meshplot!(ax1, F, V, strokewidth=0.0, color=:brown)
end

screen = display(GLMakie.Screen(), fig)     