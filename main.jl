module JuLattice

############################
## Main file for JuLattice #
############################
include("src/Logger.jl")
include("src/ConfigReader.jl")
include("src/Simulation_SetUp.jl")
include("src/Object_SetUp.jl")
include("src/Boundary_Conditions_SetUp.jl")
include("src/Plotter.jl")

using Revise, GLMakie
using .Plotter, .Logger, .ConfigReader, .Simulation_SetUp, .Object_SetUp, .Boundary_Conditions_SetUp #, .SimulationCore

## User Settings
# Domain Settings
length_X = 4;              # m
length_Y = 1;              # m 

# Object definition
Object = "Cylinder";
Radius   = 0.1    # m
Position = [1, 0.5] # m

# Fluid Settings 
Fluid_Density = 1000.0;       # kg/m^3
Inflow_Velocity = 0.4;      # m/s
Kinematic_Viscosity = 0.001; # m^2/s 

# Simulation Settings
Simulation_Time = 8000;     # s
delta_x = 0.01;             # discretisation in time and space
τ = 0.65;

# Boundary conditions
Left_BC     = "Velocity"
Left_BC_Velocity = Inflow_Velocity

Right_BC    = "ZeroGradient"
Top_BC      = "NoSlip"
Bottom_BC   = "NoSlip"

# Plot Requests
Plotvx = true;
Plotvy = true;
Plotvorticity = true;

if !@isdefined Left_BC_Velocity
    Left_BC_Velocity = 0
    if Left_BC == "Velocity"
    warning("Left BC defined as 'Velocity' but no value provided.")
    end
end
if !@isdefined Right_BC_Velocity
    Right_BC_Velocity = 0
    if Right_BC == "Velocity"
    warning("Right BC defined as 'Velocity' but no value provided.")
    end 
end
if !@isdefined Top_BC_Velocity
    Top_BC_Velocity = 0
    if Top_BC == "Velocity"
    warning("Top BC defined as 'Velocity' but no value provided.")
    end
end
if !@isdefined Bottom_BC_Velocity
    Bottom_BC_Velocity = 0
    if Bottom_BC == "Velocity"
    warning("Bottom BC defined as 'Velocity' but no value provided.")
    end
end

### Verify User Input ###
Verification(length_X, length_Y, Object, Radius, Position, Fluid_Density, Inflow_Velocity, 
Kinematic_Viscosity, Simulation_Time, delta_x, τ, 
Left_BC, Right_BC, Top_BC, Bottom_BC, Left_BC_Velocity, Right_BC_Velocity, Top_BC_Velocity, Bottom_BC_Velocity, 
Plotvx, Plotvy, Plotvorticity)

config = ReadConfig(length_X, length_Y, Object, Radius, Position, Fluid_Density, Inflow_Velocity, 
                    Kinematic_Viscosity, Simulation_Time, delta_x, τ, 
                    Left_BC, Right_BC, Top_BC, Bottom_BC, Left_BC_Velocity, Right_BC_Velocity, Top_BC_Velocity, Bottom_BC_Velocity, 
                    Plotvx, Plotvy, Plotvorticity)



#### Run Simulation #####
Log_Simulation_Header()

## Setup and initialization
simulation, fluid   = Set_Simulation_Params(config)

grid                = Create_Grid(simulation)
mutable_grid        = Initialize_Distributions(simulation, fluid, grid)

object              = Create_Object(config, simulation, fluid, grid)

boundary_conditions = Create_Boundary_Conditions(config, fluid, object)

#Log 
Log_Discretization_Settings(simulation.delta_x, simulation.delta_t, object.lattice_Reynolds, grid.gridX)
#=
## Plot Stuff to be moved 
if any((Plotvorticity, Plotvx, Plotvy))
    if Plotvorticity==true 
        vorticity, vorticity_obs, text_obj, step_text, fig_vorticity = Create_Plot(gridlengthX, gridlengthY)
        screen1 = GLMakie.Screen()
        GLMakie.display(screen1, fig_vorticity)
    end
    if Plotvx==true 
        velocityX_obs, text_obj_vx, step_text_vx, fig_vx = Create_Plot(gridlengthX, gridlengthY, velocityX, "X")
        screen2 = GLMakie.Screen(; position = (600, 0))
        GLMakie.display(screen2, fig_vx)
    end
    if Plotvy==true 
        velocityY_obs, text_obj_vy, step_text_vy, fig_vy = Create_Plot(gridlengthX, gridlengthY, velocityY, "Y")
        screen3 = GLMakie.Screen()
        display(screen3, fig_vy)
    end
end
##############################
println("#################################")
println("Starting Simulation:")

Run_Simulation!(simulation, fluid, grid, mutable_grid, object, boundary_conditions)

=#

Log_Simulation_Tail()

end