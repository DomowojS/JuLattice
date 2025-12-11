############################
## Main file for JuLattice #
############################
include("src/Plotter.jl")
include("src/Logger.jl")

using MeshGrid, GLMakie
using .Plotter, .Logger

##-------- User Settings --------##
# Domain Settings
length_X = 2              # m
length_Y = 0.5            # m 
length_Z = 1              # m

# sphere Definition
Radius   = 0.1    # m
Position = [0.5, 0.25, 0.5] # m [x,y,z]

# Fluid Settings 
Fluid_Density = 1000.0;       # kg/m^3
Inflow_Velocity = 0.4;      # m/s
Kinematic_Viscosity = 0.001; # m^2/s 

# Simulation Settings
# Simulation_Time = 8000;     # s
Simulation_Time = 5;     # s
delta_x = 0.01;             # discretisation in time and space
τ = 0.55 #0.65;
#Re = (Inflow_Velocity .* Radius)/Kinematic_Viscosity;
Re = (Inflow_Velocity .* 2 .* Radius)/Kinematic_Viscosity;
Re_Log=floor(Int,Re)

# Plot Requests (Flags)
Plotvx = false;
Plotvy = false;
Plotvz = false;
Plotvorticity = true;


####-------- Run Simulation --------#####
Log_Simulation_Header()

##-------- Compute timestep from relaxation time --------##
# Time step from relaxation time
lattice_speedOfSound = 1 / √3; #bleibt gleich bei 3D
delta_t = ((τ - 0.5) * lattice_speedOfSound^2 * delta_x^2) / Kinematic_Viscosity

##-------- Convert user settings to lattice units --------##
# Domain
gridlengthX = ceil(Int, length_X / delta_x);
gridlengthY = ceil(Int, length_Y / delta_x);
gridlengthZ = ceil(Int, length_Z / delta_x);

# Sphere
sphere_radius  = Radius/delta_x;
sphere_position = Position ./ delta_x;

# Fluid
fluiddensity = Fluid_Density
#fluiddensity = 100;
lattice_inflow_velocity = Inflow_Velocity * (delta_t / delta_x);    # m/s  * s/m -> [-]
lattice_viscosity = lattice_speedOfSound^2 * (τ -0.5);              #(m/s)^2 * s = m^2/s

# ReynoldsCheck
lattice_Re = (lattice_inflow_velocity .* sphere_radius)/lattice_viscosity; #Re_lattice = U*R/v -> sollte Re entsprechen weil Größen skaliert wurden
lattice_Re_Log=floor(Int,lattice_Re)

# Log 
Log_Discretization_Settings(delta_x, delta_t, lattice_Re_Log)

# Simulation Settings
simulationTime = ceil(Int, Simulation_Time / delta_t);

Q = 19; #D3Q19
#link Hand Koordinatensystem
# velocity_vector = [
#     [0, 1, -1, 0,  0, 0,  0, 1, -1,  1, -1, 0,  0,  0,  0, 1, -1, -1,  1],        #x
#     [0, 0,  0, 0,  0, 1, -1, 0,  0,  0,  0, 1, -1, -1,  1, 1, -1,  1, -1],        #y
#     [0, 0,  0, 1, -1, 0,  0, 1, -1, -1,  1, 1, -1,  1, -1, 0,  0,  0,  0]         #z
# ];

#rechte Hand Koordinatensystem
velocity_vector = [
    [0, 1, -1, 0,  0,  0, 0, 1, -1,  1, -1,  0, 0, 0,  0,  1, -1, -1, 1],        #x
    [0, 0,  0, 0,  0, -1, 1, 0,  0,  0,  0, -1, 1, 1, -1, -1,  1, -1, 1],        #y
    [0, 0,  0, 1, -1,  0, 0, 1, -1, -1,  1,  1,-1, 1, -1,  0,  0,  0, 0]         #z
];


#define velocity vectors for all directions
velocity_vector_x = reshape(velocity_vector[1,],1,1,1,Q)
velocity_vector_y = reshape(velocity_vector[2,],1,1,1,Q)
velocity_vector_z = reshape(velocity_vector[3,],1,1,1,Q)


weights = [3/9,                                      #ruhe (1)
            1/18, 1/18, 1/18, 1/18, 1/18, 1/18,      #axial (6)
            1/36, 1/36, 1/36, 1/36, 1/36, 1/36,      #diagonal (12)
            1/36, 1/36, 1/36, 1/36, 1/36, 1/36];     #-------------19

weights = reshape(weights, 1, 1, 1, Q);

# create grid
gridX, gridY, gridZ = meshgrid(1:gridlengthX, 1:gridlengthY, 1:gridlengthZ);

#Swap of Y and X axis: (Y,X,Z) -> (X,Y,Z)
if size(gridX) == (gridlengthY, gridlengthX, gridlengthZ) #check for format of grids
    gridX = permutedims(gridX, (2,1,3))
    gridY = permutedims(gridY, (2,1,3))
    gridZ = permutedims(gridZ, (2,1,3))
end    

midY = ceil(Int, gridlengthY/2)
midZ = ceil(Int, gridlengthZ/2)

# #---- Debug
# @show gridlengthX, gridlengthY, gridlengthZ
# @show size(gridX), size(gridY), size(gridZ)
# @show gridX[1,1,1], gridX[end,1,1]
# @show gridY[1,1,1], gridY[1,end,1]
# @show gridZ[1,1,1], gridZ[1,1,end]
# flush(stdout)
# error("DEBUG STOP")
# #---- Debug

# create object indetifier
sphere = (gridX .- sphere_position[1]).^2 + (gridY .- sphere_position[2]).^2 + (gridZ .- sphere_position[3]).^2 .< sphere_radius.^2

# create boundary indetifiers
walls = gridY .== 1 .| gridY .== gridlengthY .| gridZ .==1 .| gridZ .== gridlengthZ;
inlet = gridX .== 1;
outlet = gridX .== gridlengthX;

# Initialize distributions arrays 3D
distributions = ones(gridlengthX, gridlengthY, gridlengthZ, Q) .+ 0.01*rand(gridlengthX, gridlengthY, gridlengthZ, Q);
distributions[:,:,:,2 ] .+= 2 .* (1 .+ 0.2 .* cos.(2 .* π .* gridX ./ gridlengthX .*4));
distributions_equilibrium = ones(gridlengthX, gridlengthY, gridlengthZ, Q);

# Initialize macroscopic density and scale distribution 3D
densityGrid = sum(distributions, dims=4);
distributions .*= fluiddensity ./ densityGrid;

# Initialize macroscopic velocity arrays
velocityX   = zeros(gridlengthX, gridlengthY, gridlengthZ);
velocityY   = zeros(gridlengthX, gridlengthY, gridlengthZ);
velocityZ   = zeros(gridlengthX, gridlengthY, gridlengthZ);

#Initialise macroscopic vorticity arrays
omegaX = zeros(gridlengthX, gridlengthY, gridlengthZ)
omegaY = zeros(gridlengthX, gridlengthY, gridlengthZ)
omegaZ = zeros(gridlengthX, gridlengthY, gridlengthZ)
omegaMag = zeros(gridlengthX, gridlengthY, gridlengthZ)
# Initialise dotproduct array 
dotprod_velocities = zeros(gridlengthX, gridlengthY, gridlengthZ, Q);

#Plot calls
if any((Plotvorticity, Plotvx, Plotvy, Plotvz))
    # #3D Plot
    # if Plotvorticity == true
    #     omegaMag_obs, step_text_omega, fig_omega = Create_Plot3D(gridlengthX, gridlengthY, gridlengthZ, omegaMag; title="|ω|")
    #     screenOmega = GLMakie.Screen()
    #     display(screenOmega, fig_omega)
    # end
    
    # Bildschirmgröße ermitteln
    #Default
    screen_width = 1920
    screen_height = 1080

    try
        monitor = GLMakie.GLFW.GetPrimaryMonitor()
        mode = GLMakie.GLFW.GetVideoMode(monitor)
        global screen_width = mode.width
        global screen_height = mode.height
    catch
    end

    gap = 20
    window_width = Int(floor((screen_width - 3*gap) / 2))
    window_height = Int(floor((screen_height - 3*gap) / 2))

    ###DEBUG###
    println("Screen Resolution: $(screen_width) x $(screen_height)")

    if Plotvx==true
        #xy slice at z=midZ
        vx_xy_obs, step_text_vx_xy, fig_vx_xy = Create_Plot_XY(gridlengthX, gridlengthY, velocityX[:,:,midZ]; width=window_width, height=window_height, title="v_x at z=$(midZ)")
        screen_vx_xy = GLMakie.Screen(position=(gap, gap))
        display(screen_vx_xy, fig_vx_xy)

        #xz slice at y=midY
        vx_xz_obs, step_text_vx_xz, fig_vx_xz = Create_Plot_XZ(gridlengthX, gridlengthZ, velocityX[:,midY,:]; width=window_width, height=window_height, title="v_x at y=$(midY)")
        screen_vx_xz = GLMakie.Screen(position=(window_width + 2*gap, gap))
        display(screen_vx_xz, fig_vx_xz)
    end

    if Plotvorticity == true
        #xy slice at z=midZ
        omega_xy_obs, step_text_omega_xy, fig_omega_xy = Create_Vorticity_XY(gridlengthX, gridlengthY, omegaMag[:,:,midZ]; width=window_width, height=window_height, title = "|ω| at z=$(midZ)")
        screen_omega_xy = GLMakie.Screen(position=(gap, window_height + 2*gap))
        display(screen_omega_xy, fig_omega_xy)

        #xz slice at y=midY
        omega_xz_obs, step_text_omega_xz, fig_omega_xz = Create_Vorticity_XZ(gridlengthX, gridlengthZ, omegaMag[:,midY,:]; width=window_width, height=window_height, title = "|ω| at y=$(midY)")
        screen_omega_xz = GLMakie.Screen(position=(window_width + 2*gap, window_height + 2*gap))
        display(screen_omega_xz, fig_omega_xz)    
    end


    if Plotvy==true 
        velocityY_obs, text_obj_vy, step_text_vy, fig_vy = Create_Plot(gridlengthX, gridlengthY, velocityY, "Y")
        screen3 = GLMakie.Screen()
        display(screen3, fig_vy)
    end

    if Plotvz==true
        velocityZ_obs, text_obj_vz, step_text_vz, fig_vz = Create_Plot(gridlengthX, gridlengthZ, velocityZ, "Z")
        screen4 = GLMakie.Screen()
        display(screen4, fig_vz)
    end

end

println("#################################")
println("Starting Simulation:")
# Run Simulation Loop
for i in 1:simulationTime

    # Get Macroscopic values 3D
    global densityGrid = sum(distributions, dims=4);
    velocityX .= (1 ./ densityGrid) .* sum(distributions.*velocity_vector_x, dims=4); 
    velocityY .= (1 ./ densityGrid) .* sum(distributions.*velocity_vector_y, dims=4); 
    velocityZ .= (1 ./ densityGrid) .* sum(distributions.*velocity_vector_z, dims=4); 

    #DEBUGGING############
    if i % 20 == 0
        @views begin
            xy_min, xy_max = extrema(velocityX[:,:,midZ])
            xz_min, xz_max = extrema(velocityX[:,midY,:])
        end

    println("vx XY min/max : ($(xy_min), $(xy_max)) | XZ min/max: ($(xz_min), $(xz_max)) at step $i")
    end
    #DEBUGGING############


    ## Apply Collision 3D
    # Compute equilibrium state
    dotprod_velocities .= (velocity_vector_x .* velocityX) .+ (velocity_vector_y .* velocityY) .+ (velocity_vector_z .* velocityZ);
    distributions_equilibrium .= weights .* densityGrid .*(1 .+ 3 .*dotprod_velocities .+ 4.5 .*dotprod_velocities.^2 .- 1.5 .*(velocityX.^2 .+ velocityY.^2 + velocityZ.^2));
    # Relax towards equilibrium
    distributions .+= -(1/τ) .* (distributions .- distributions_equilibrium);

    # Stream 
    for j in 1:Q
        distributions[:,:,:,j] = circshift(distributions[:,:,:,j], (velocity_vector_x[j], velocity_vector_y[j], velocity_vector_z[j]))
    end

    # ------------------------ vorerst ohne ------------------------ 
    # ## Apply Boundary conditions
    # #Inlet velocity bc (unknown: f_1, f_8, f_9)
    # densityGrid[inlet, :] .= (sum(distributions[inlet, [1,3,5]], dims=2).+ 2 .*sum(distributions[inlet, [2,6,7]], dims=2)) ./ (1-lattice_inflow_velocity)
    # distributions[inlet, 4] .= distributions[inlet, 2] .+ (2/3 .* densityGrid[inlet,:] .* lattice_inflow_velocity)
    # distributions[inlet, 8] .= distributions[inlet, 6] .+ (1/6 .* densityGrid[inlet,:] .* lattice_inflow_velocity) .- (1/2 .* (distributions[inlet, 3] .- distributions[inlet, 5]))
    # distributions[inlet, 9] .= distributions[inlet, 7] .+ (1/6 .* densityGrid[inlet,:] .* lattice_inflow_velocity) .+ (1/2 .* (distributions[inlet, 3] .- distributions[inlet, 5]))
    

    # #Outlet zero gradient bc
    # distributions[outlet, [4, 8, 9]] .= distributions[gridlengthX-1, :, [4, 8, 9]]
    # ------------------------ vorerst ohne ------------------------ 

    #No Slip Walls 3D
    distributions[walls, 1:Q] .= distributions[walls, [1,3,2,5,4,7,6,9,8,11,10,13,12,15,14,17,16,19,18]];
    
    # Apply object boundary condition 3D
    distributions[sphere, 1:Q] .= distributions[sphere, [1,3,2,5,4,7,6,9,8,11,10,13,12,15,14,17,16,19,18]];

        # Plot of the field
        if ((i % 10 == 0)) || (i == simulationTime)
            # Set velocities inside the sphere to zero
            velocityX[sphere] .= NaN
            velocityY[sphere] .= NaN
            velocityZ[sphere] .= NaN

            # # Compute vorticity 3D
            #vorticity
            dv_dx = (circshift(velocityY, (-1,0,0)) .- circshift(velocityY, (1,0,0))) ./ 2
            dw_dx = (circshift(velocityZ, (-1,0,0)) .- circshift(velocityZ, (1,0,0))) ./ 2

            du_dy = (circshift(velocityX, (0,-1,0)) .- circshift(velocityX, (0,1,0))) ./ 2
            dw_dy = (circshift(velocityZ, (0,-1,0)) .- circshift(velocityZ, (0,1,0))) ./ 2
            
            du_dz = (circshift(velocityX, (0,0,-1)) .- circshift(velocityX, (0,0,1))) ./ 2
            dv_dz = (circshift(velocityY, (0,0,-1)) .- circshift(velocityY, (0,0,1))) ./ 2

            omegaX .= dw_dy .- dv_dz
            omegaY .= du_dz .- dw_dx
            omegaZ .= dv_dx .-du_dy
            omegaMag .= sqrt.(omegaX.^2 .+ omegaY.^2 .+ omegaZ.^2)
            
            if ((i % 100 == 0)) || (i == simulationTime)
                Log_Simulation_Runtime(i, simulationTime)
            end
            # Update the observables
            if Plotvorticity == true   
                                   
                #Mask inlet outlet
                omegaMag[inlet] .= 0.0
                omegaMag[outlet] .= 0.0
                # Mask the cylinder region
                omegaMag[sphere] .= NaN
                
                # #3D update: 
                # omegaMag_obs[] = Float32.(omegaMag)
                # step_text_omega[] = "Time step: $i, $(floor(Int, i*delta_t))s"
                
                #2D update:
                #xy at midZ
                omega_xy_obs[] = copy(omegaMag[:,:,midZ])
                step_text_omega_xy[] = "Time step: $i, $(floor(Int, i*delta_t))s" 
                #xz at midY
                omega_xz_obs[] = copy(omegaMag[:,midY,:])
                step_text_omega_xz[] = "Time step: $i, $(floor(Int, i*delta_t))s"

            end

            if Plotvx==true
                #refresh xy slice at z=midZ
                vx_xy_obs[] = copy(velocityX[:,:,midZ])
                step_text_vx_xy[] = "Time step:$i, $(floor(Int, i*delta_t))s"

                #refresh xz slice at y=midY
                vx_xz_obs[] = copy(velocityX[:,midY,:])
                step_text_vx_xz[] = "Time step:$i, $(floor(Int, i*delta_t))s"
            end

            if Plotvy==true 
                velocityY_obs[] = copy(velocityY) 
                step_text_vy[] = "Time step: $i, $(floor(Int, i*delta_t))s"
            end

            yield()
            sleep(0.05)
        end

end

Log_Simulation_Tail()