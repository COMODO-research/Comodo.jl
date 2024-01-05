#Visualize mesh
fig = Figure()

stepRange = 0:1:5
sl_step = Slider(fig[2, 1], range = stepRange, startvalue = 0)

y = lift(sl_step.value) do stepIndex
    return stepIndex
end

titleString = lift(sl_step.value) do stepIndex
  "Step: "*string(stepIndex)
end

ax=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", 
title = titleString,limits = (-2, 2, 0, 5, -2, 2))
hp=scatter!(0,y,markersize = 25)

g = Observable(false)

# hb=Button(fig, label = "Start")
# on(hb.clicks) do n
obs_func = on(g) do _
    @async for s ∈ stepRange
        set_close_to!(sl_step, s)      
        sleep(0.1)
    end
end    
# end
ht = Toggle(fig, active = false)

connect!(g, ht.active)

fig
