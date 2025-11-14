using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra

GLMakie.closeall()

nSteps = 50
xr,yr,zr = ntuple(_->range(0,4*pi,nSteps),3)
level = 0.0
cap = true

# Visualization
fig = Figure(size=(800,800))

typeSet = (:P,   # Schwarz P-surface
           :D,   # Schwarz D-surface        
           :D2,  # Schwarz D-surface        
           :N,   # Neovius
           :G,   # Gyroid 
           :IWP, # Schoen I-graph-wrapped package 
           :FRD, # Schoen faces of rhombic-dodecahedron
           :S,   # Fish-Kocher S-sheet
           :HG)  # Lidinoid

s = ceil(Int,sqrt(length(typeSet))) # Visualisation grid size         
for (i,triplyPeriodicType) in enumerate(typeSet)

    A = [triplyperiodicminimal((x,y,z), triplyPeriodicType) for x in xr, y in yr, z in zr]    
    F1, V1 = getisosurface(A; x = xr, y = yr, z = zr, level = level, cap = cap, padValue=1e8)      

    indPlot = CartesianIndices((s,s))[i]
    ax1 = AxisGeom(fig[indPlot[1], indPlot[2]]; title="type=:$triplyPeriodicType", limits=(minimum(xr),maximum(xr),minimum(yr),maximum(yr),minimum(zr),maximum(zr)))
    hp1 = meshplot!(ax1, F1, V1, color=:brown, strokewidth=0.0)
end
screen = display(GLMakie.Screen(), fig)