using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra

GLMakie.closeall()

# Visualization
fig = Figure(size=(1400,1000))
cartInd = CartesianIndices((3,4))
for (i_plot, tpms_type) in enumerate( (:P, :D, :D2, :N, :G, :IWP, :FRD, :S, :HG) )       
    level = 0.0
    F1, V1 = tpms(tpms_type; x=range(0, 4*pi, 50), level=-level, cap=true, side=:positive)        
    F2, V2 = tpms(tpms_type; x=range(0, 4*pi, 50), level=level, cap=true, side=:negative)        

    ij = cartInd[i_plot]
    ax1 = AxisGeom(fig[ij[1], ij[2]]; title=string(tpms_type))
    hp1 = meshplot!(ax1, F1, V1, strokewidth=0.0, color=:brown)
    hp2 = meshplot!(ax1, F2, V2, strokewidth=0.0, color=(:white,0.25), transparency=true)
end

screen = display(GLMakie.Screen(), fig)     