using XML
using Printf

######

V = [0.0 0.0 0.0; 
     1.0 0.0 0.0;
     1.0 1.0 0.0;
     0.0 1.0 0.0;
     0.0 0.0 1.0; 
     1.0 0.0 1.0;
     1.0 1.0 1.0;
     0.0 1.0 1.0;]

E =[1 2 3 4 5 6 7 8]

function aen(main_node::Node,sub_node_name::String, args...; kwargs...)
    if isnothing(args)
        sub_node = XML.Element(sub_node_name; kwargs...)    
    else
        sub_node = XML.Element(sub_node_name, args...; kwargs...)    
    end
    push!(main_node, sub_node)
    return sub_node
end

function febIni()
    # Define document level
    doc = XML.Document()

    # Define declaration: <?xml version="1.0" encoding="UTF-8"?>
    declarationNode = XML.Declaration(version = "1.0", encoding = "UTF-8")
    push!(doc, declarationNode)

    # Define febio_spec node: <febio_spec version="4.0">
    febio_spec_node = XML.Element("febio_spec"; version = "4.0")
    push!(doc, febio_spec_node)

    # Add comment: <!--Created using Julia-->
    comment_node = XML.Comment("Created using Julia")
    push!(febio_spec_node, comment_node)
    
    return doc,febio_spec_node
end


######

# Define file name
filenameOut = "febioInputFile_01.feb"

doc,febio_spec_node = febIni()

aen(febio_spec_node,"Module"; type = "solid") # Define Module node: <Module type="solid"/>

control_node = aen(febio_spec_node,"Control") # Define Control node: <Control>
    aen(control_node,"analysis","STATIC")               
    aen(control_node,"time_steps",40)
    aen(control_node,"step_size",1e-2)
    aen(control_node,"plot_zero_state",1)
    aen(control_node,"plot_range",@sprintf("%.2f, %.2f",0,-1))
    aen(control_node,"plot_level","PLOT_MAJOR_ITRS")
    aen(control_node,"plot_stride",1)
    aen(control_node,"output_level","OUTPUT_MAJOR_ITRS")
    aen(control_node,"adaptor_re_solve",1)

time_stepper_node = aen(control_node,"time_stepper"; type = "default")
    aen(time_stepper_node,"max_retries",5)
    aen(time_stepper_node,"opt_iter",35)
    aen(time_stepper_node,"dtmin",2.5e-4)
    aen(time_stepper_node,"dtmax",5e-2)
    aen(time_stepper_node,"aggressiveness",0)
    aen(time_stepper_node,"cutback",5e-1)
    aen(time_stepper_node,"dtforce",0)

solver_node = aen(control_node,"solver"; type = "solid")
    aen(solver_node,"symmetric_stiffness",1)
    aen(solver_node,"equation_scheme",1)
    aen(solver_node,"equation_order","default")
    aen(solver_node,"optimize_bw",0)
    aen(solver_node,"lstol",9e-1)
    aen(solver_node,"lsmin",1e-2)
    aen(solver_node,"lsiter",5)
    aen(solver_node,"max_refs",70)
    aen(solver_node,"check_zero_diagonal",0)
    aen(solver_node,"zero_diagonal_tol",0)
    aen(solver_node,"force_partition",0)
    aen(solver_node,"reform_each_time_step",1)
    aen(solver_node,"reform_augment",0)
    aen(solver_node,"diverge_reform",1)
    aen(solver_node,"min_residual",1e-20)
    aen(solver_node,"max_residual",0)
    aen(solver_node,"dtol",1e-3)
    aen(solver_node,"etol",1e-2)
    aen(solver_node,"rtol",0)
    aen(solver_node,"rhoi",0)
    aen(solver_node,"alpha",1)
    aen(solver_node,"beta",2.5e-01)
    aen(solver_node,"gamma",5e-01)
    aen(solver_node,"logSolve",0)
    aen(solver_node,"arc_length",0)
    aen(solver_node,"arc_length_scale",0)
qn_method_node = aen(solver_node,"qn_method"; type = "BFGS")
    aen(qn_method_node,"max_ups",0)
    aen(qn_method_node,"max_buffer_size",0)
    aen(qn_method_node,"cycle_buffer",0)
    aen(qn_method_node,"cmax",0)

Globals_node   = aen(febio_spec_node,"Globals")

Constants_node = aen(Globals_node,"Constants")
    aen(Constants_node,"R",8.3140000e-06)
    aen(Constants_node,"T",298)
    aen(Constants_node,"F",9.6485000e-05)

Material_node = aen(febio_spec_node,"Material")

material_node = aen(Material_node,"material"; id = 1, name="Material1", type="neo-Hookean")
    aen(material_node,"E",1)
    aen(material_node,"v",0.4)

Mesh_node = aen(febio_spec_node,"Mesh")

Nodes_node = aen(Mesh_node,"Nodes"; name="nodeSet_all")
    for q ∈ 1:size(V,1)
        aen(Nodes_node,"node",@sprintf("%.2f, %.2f, %.2f",V[q,1],V[q,2],V[q,3]); id = q)
    end
    
# Elements
Elements_node = aen(Mesh_node,"Elements"; name="Part1", type="hex8")
    for q ∈ 1:size(E,1)
        aen(Elements_node,"elem",@sprintf("%i, %i, %i, %i, %i, %i, %i, %i",E[q,1],E[q,2],E[q,3],E[q,4],E[q,5],E[q,6],E[q,7],E[q,8]); id = q)
    end
# Node sets
aen(Mesh_node,"NodeSet",@sprintf("%i, %i, %i, %i",E[1,1],E[1,2],E[1,3],E[1,4]); name="bcSupportList")
aen(Mesh_node,"NodeSet",@sprintf("%i, %i, %i, %i",E[1,5],E[1,6],E[1,7],E[1,8]); name="bcPrescribeList")

MeshDomains_node = aen(febio_spec_node, "MeshDomains")
    aen(MeshDomains_node,"SolidDomain"; mat = "Material1", name="Part1")

Boundary_node = aen(febio_spec_node, "Boundary")

bc_node = aen(Boundary_node,"bc"; name="zero_displacement_xyz", node_set="bcSupportList", type="zero displacement")
    aen(bc_node,"x_dof",1)
    aen(bc_node,"y_dof",1)
    aen(bc_node,"z_dof",1)

bc_node2 = aen(Boundary_node,"bc"; name="prescribed_disp_x", node_set="bcPrescribeList", type="prescribed displacement")
    aen(bc_node2,"dof","x")
    aen(bc_node2,"value",1; lc=@sprintf("%i",1))
    aen(bc_node2,"relative",@sprintf("%i",0))

bc_node3 = aen(Boundary_node,"bc"; name="prescribed_disp_y", node_set="bcPrescribeList", type="prescribed displacement")
    aen(bc_node3,"dof","y")
    aen(bc_node3,"value",0; lc=@sprintf("%i",1))
    aen(bc_node3,"relative",@sprintf("%i",0))

bc_node4 = aen(Boundary_node,"bc"; name="prescribed_disp_z", node_set="bcPrescribeList", type="prescribed displacement")
    aen(bc_node4,"dof","z")
    aen(bc_node4,"value",0; lc=@sprintf("%i",1))
    aen(bc_node4,"relative",@sprintf("%i",0))

LoadData_node = aen(febio_spec_node,"LoadData")

load_controller_node = aen(LoadData_node,"load_controller"; id=1, name="LC_1", type="loadcurve")
    aen(load_controller_node,"interpolate","LINEAR")
    
points_node = aen(load_controller_node,"points")
    aen(points_node,"pt",@sprintf("%.2f, %.2f",0,0))
    aen(points_node,"pt",@sprintf("%.2f, %.2f",1,1))

Output_node = aen(febio_spec_node,"Output")

plotfile_node = aen(Output_node,"plotfile"; type="febio")
    aen(plotfile_node,"var"; type="displacement")
    aen(plotfile_node,"var"; type="stress")
    aen(plotfile_node,"var"; type="relative volume")
    aen(plotfile_node,"var"; type="reaction forces")
    aen(plotfile_node,"var"; type="contact pressure")
    aen(plotfile_node,"compression",@sprintf("%i",0))

#######

# Write to XML file
XML.write(filenameOut, doc)

FEBioPath = "/home/kevin/FEBioStudio/bin/febio4"
run_filename = filenameOut
runCommand = `nice "$FEBioPath" "$run_filename"`
run(runCommand);