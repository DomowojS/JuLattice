############################
## Main file for JuLattice #
############################
include("src/Plotter.jl")
include("src/Logger.jl")

using MeshGrid, GLMakie
using .Plotter, .Logger

####################################  Initialize  ####################################
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
τ = 0.65 #0.65;
#Re = (Inflow_Velocity .* Radius)/Kinematic_Viscosity;
Re = (Inflow_Velocity .* 2 .* Radius)/Kinematic_Viscosity;
Re_Log=floor(Int,Re)

# Plot Requests (Flags)
Plotvx = true;
Plotvy = false;
Plotvz = false;
Plotvorticity = false;


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
# D3Q19
# f000 = rest (0,0,0)
# fm00, fp00 = x-axis (±1,0,0)
# f0m0, f0p0 = y-axis (0,±1,0)
# f00m, f00p = z-axis (0,0,±1)
# fmm0, fmp0, fpm0, fpp0 = xy-plane edges
# fm0m, fm0p, fp0m, fp0p = xz-plane edges
# f0mm, f0mp, f0pm, f0pp = yz-plane edges

#Define arrays for each direction D3Q9
# f000 = rest (0,0,0)
f000 = zeros(gridlengthX, gridlengthY, gridlengthZ)
# fm00, fp00 = x-axis (±1,0,0)
fm00 = zeros(gridlengthX, gridlengthY, gridlengthZ)
fp00 = zeros(gridlengthX, gridlengthY, gridlengthZ)
# f0m0, f0p0 = y-axis (0,±1,0)
f0m0 = zeros(gridlengthX, gridlengthY, gridlengthZ)
f0p0 = zeros(gridlengthX, gridlengthY, gridlengthZ)
# f00m, f00p = z-axis (0,0,±1)
f00m = zeros(gridlengthX, gridlengthY, gridlengthZ)
f00p = zeros(gridlengthX, gridlengthY, gridlengthZ)
# fmm0, fmp0, fpm0, fpp0 = xy-plane edges
fmm0 = zeros(gridlengthX, gridlengthY, gridlengthZ)
fmp0 = zeros(gridlengthX, gridlengthY, gridlengthZ)
fpm0 = zeros(gridlengthX, gridlengthY, gridlengthZ)
fpp0 = zeros(gridlengthX, gridlengthY, gridlengthZ)
# fm0m, fm0p, fp0m, fp0p = xz-plane edges
fm0m = zeros(gridlengthX, gridlengthY, gridlengthZ)
fm0p = zeros(gridlengthX, gridlengthY, gridlengthZ)
fp0m = zeros(gridlengthX, gridlengthY, gridlengthZ)
fp0p = zeros(gridlengthX, gridlengthY, gridlengthZ)
# f0mm, f0mp, f0pm, f0pp = yz-plane edges
f0mm = zeros(gridlengthX, gridlengthY, gridlengthZ)
f0mp = zeros(gridlengthX, gridlengthY, gridlengthZ)
f0pm = zeros(gridlengthX, gridlengthY, gridlengthZ)
f0pp = zeros(gridlengthX, gridlengthY, gridlengthZ)

#Define array for each direction after Collision+stream (S)
# f000 = rest (0,0,0)
f000S = zeros(gridlengthX, gridlengthY, gridlengthZ)
# fm00, fp00 = x-axis (±1,0,0)
fm00S = zeros(gridlengthX, gridlengthY, gridlengthZ)
fp00S = zeros(gridlengthX, gridlengthY, gridlengthZ)
# f0m0, f0p0 = y-axis (0,±1,0)
f0m0S = zeros(gridlengthX, gridlengthY, gridlengthZ)
f0p0S = zeros(gridlengthX, gridlengthY, gridlengthZ)
# f00m, f00p = z-axis (0,0,±1)
f00mS = zeros(gridlengthX, gridlengthY, gridlengthZ)
f00pS = zeros(gridlengthX, gridlengthY, gridlengthZ)
# fmm0, fmp0, fpm0, fpp0 = xy-plane edges
fmm0S = zeros(gridlengthX, gridlengthY, gridlengthZ)
fmp0S = zeros(gridlengthX, gridlengthY, gridlengthZ)
fpm0S = zeros(gridlengthX, gridlengthY, gridlengthZ)
fpp0S = zeros(gridlengthX, gridlengthY, gridlengthZ)
# fm0m, fm0p, fp0m, fp0p = xz-plane edges
fm0mS = zeros(gridlengthX, gridlengthY, gridlengthZ)
fm0pS = zeros(gridlengthX, gridlengthY, gridlengthZ)
fp0mS = zeros(gridlengthX, gridlengthY, gridlengthZ)
fp0pS = zeros(gridlengthX, gridlengthY, gridlengthZ)
# f0mm, f0mp, f0pm, f0pp = yz-plane edges
f0mmS = zeros(gridlengthX, gridlengthY, gridlengthZ)
f0mpS = zeros(gridlengthX, gridlengthY, gridlengthZ)
f0pmS = zeros(gridlengthX, gridlengthY, gridlengthZ)
f0ppS = zeros(gridlengthX, gridlengthY, gridlengthZ)

#Initialise macroscopic variables
rho = ones(gridlengthX, gridlengthY, gridlengthZ) .* fluiddensity
u = zeros(gridlengthX, gridlengthY, gridlengthZ)    #ux
v = zeros(gridlengthX, gridlengthY, gridlengthZ)    #uy
w = zeros(gridlengthX, gridlengthY, gridlengthZ)    #uz

#Define omega
omega = 1.0 / τ

# #create grid dimensions
# rows = gridlengthY
# cols = gridlengthX
# slices = gridlengthZ

# create grid
gridX, gridY, gridZ = meshgrid(1:gridlengthX, 1:gridlengthY, 1:gridlengthZ);

#Swap of Y and X axis: (Y,X,Z) -> (X,Y,Z)
if size(gridX) == (gridlengthY, gridlengthX, gridlengthZ) #check for format of grids
    gridX = permutedims(gridX, (2,1,3))
    gridY = permutedims(gridY, (2,1,3))
    gridZ = permutedims(gridZ, (2,1,3))
end    

# Initialise velocity arrays for plotting
velocityX = zeros(gridlengthX, gridlengthY, gridlengthZ)
velocityY = zeros(gridlengthX, gridlengthY, gridlengthZ)
velocityZ = zeros(gridlengthX, gridlengthY, gridlengthZ)

omegaX = zeros(gridlengthX, gridlengthY, gridlengthZ)
omegaY = zeros(gridlengthX, gridlengthY, gridlengthZ)
omegaZ = zeros(gridlengthX, gridlengthY, gridlengthZ)
omegaMag = zeros(gridlengthX, gridlengthY, gridlengthZ)

#Define mid-Planes for plotting
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

#Initialize distribution functions
for x in 1:gridlengthX
    for y in 1:gridlengthY
        for z in 1:gridlengthZ
            
            ux = lattice_inflow_velocity
            uy = 0.0
            uz = 0.0
            rho_init = fluiddensity
            
            # Pre-compute polynomial factors
            ux2 = ux * ux
            uy2 = uy * uy
            uz2 = uz * uz
            
            Pm_u = 1 - 3*ux + 3*ux2
            P0_u = -2 + 3*ux2
            Pp_u = 1 + 3*ux + 3*ux2
            
            Pm_v = 1 - 3*uy + 3*uy2
            P0_v = -2 + 3*uy2
            Pp_v = 1 + 3*uy + 3*uy2
            
            Pm_w = 1 - 3*uz + 3*uz2
            P0_w = -2 + 3*uz2
            Pp_w = 1 + 3*uz + 3*uz2
            
            # Push scheme: Rest particle (0,0,0) - weight 1/3
            f000[x,y,z] = rho_init * P0_u * P0_v * P0_w / 3.0
            
            # Push scheme: Face neighbors - weight 1/18
            fm00[x,y,z] = rho_init * Pm_u * P0_v * P0_w / 18.0 
            fp00[x,y,z] = rho_init * Pp_u * P0_v * P0_w / 18.0 
            
            f0m0[x,y,z] = rho_init * P0_u * Pm_v * P0_w / 18.0 
            f0p0[x,y,z] = rho_init * P0_u * Pp_v * P0_w / 18.0
            
            f00m[x,y,z] = rho_init * P0_u * P0_v * Pm_w / 18.0 
            f00p[x,y,z] = rho_init * P0_u * P0_v * Pp_w / 18.0 
            
            # Push scheme: Edge neighbors - weight 1/36
            # XY-plane edges
            fmm0[x,y,z] = rho_init * Pm_u * Pm_v * P0_w / 36.0 
            fmp0[x,y,z] = rho_init * Pm_u * Pp_v * P0_w / 36.0 
            fpm0[x,y,z] = rho_init * Pp_u * Pm_v * P0_w / 36.0 
            fpp0[x,y,z] = rho_init * Pp_u * Pp_v * P0_w / 36.0 
            
            # XZ-plane edges
            fm0m[x,y,z] = rho_init * Pm_u * P0_v * Pm_w / 36.0 
            fm0p[x,y,z] = rho_init * Pm_u * P0_v * Pp_w / 36.0
            fp0m[x,y,z] = rho_init * Pp_u * P0_v * Pm_w / 36.0 
            fp0p[x,y,z] = rho_init * Pp_u * P0_v * Pp_w / 36.0 
            
            # YZ-plane edges
            f0mm[x,y,z] = rho_init * P0_u * Pm_v * Pm_w / 36.0 
            f0mp[x,y,z] = rho_init * P0_u * Pm_v * Pp_w / 36.0 
            f0pm[x,y,z] = rho_init * P0_u * Pp_v * Pm_w / 36.0 
            f0pp[x,y,z] = rho_init * P0_u * Pp_v * Pp_w / 36.0
        end
    end
end

#Initialise fS-Arrays for the first time as f-Arrays
f000S .= f000 
fm00S .= fm00
fp00S .= fp00
f0m0S .= f0m0 
f0p0S .= f0p0 
f00mS .= f00m 
f00pS .= f00p 
fmm0S .= fmm0 
fmp0S .= fmp0 
fpm0S .= fpm0 
fpp0S .= fpp0 
fm0mS .= fm0m 
fm0pS .= fm0p 
fp0mS .= fp0m 
fp0pS .= fp0p 
f0mmS .= f0mm 
f0mpS .= f0mp 
f0pmS .= f0pm
f0ppS .= f0pp

#Density check vielleicht?
# for x in 1:gridlengthX
#     for y in 1:gridlengthY
        
#         rho_check = f00[x,y] + fm0[x,y] + fp0[x,y] + f0m[x,y] + f0p[x,y] +
#                     fmm[x,y] + fmp[x,y] + fpm[x,y] + fpp[x,y]
        
#         rho_fluiddensity = fluiddensity

#         if abs(rho_check - rho_fluiddensity) > 0.01 * rho_fluiddensity
#             println("=== DENSITY CHECK ===")
#             println("rho_check: $rho_check , fluiddensity: $rho_fluiddensity")
#         end
#     end
# end


# # Initialize distributions arrays 3D
# distributions = ones(gridlengthX, gridlengthY, gridlengthZ, Q) .+ 0.01*rand(gridlengthX, gridlengthY, gridlengthZ, Q);
# distributions[:,:,:,2 ] .+= 2 .* (1 .+ 0.2 .* cos.(2 .* π .* gridX ./ gridlengthX .*4));
# distributions_equilibrium = ones(gridlengthX, gridlengthY, gridlengthZ, Q);

# # Initialize macroscopic density and scale distribution 3D
# densityGrid = sum(distributions, dims=4);
# distributions .*= fluiddensity ./ densityGrid;

# # Initialize macroscopic velocity arrays
# velocityX   = zeros(gridlengthX, gridlengthY, gridlengthZ);
# velocityY   = zeros(gridlengthX, gridlengthY, gridlengthZ);
# velocityZ   = zeros(gridlengthX, gridlengthY, gridlengthZ);

# #Initialise macroscopic vorticity arrays
# omegaX = zeros(gridlengthX, gridlengthY, gridlengthZ)
# omegaY = zeros(gridlengthX, gridlengthY, gridlengthZ)
# omegaZ = zeros(gridlengthX, gridlengthY, gridlengthZ)
# omegaMag = zeros(gridlengthX, gridlengthY, gridlengthZ)
# # Initialise dotproduct array 
# dotprod_velocities = zeros(gridlengthX, gridlengthY, gridlengthZ, Q);

#Plot calls
if any((Plotvorticity, Plotvx, Plotvy, Plotvz))
    # #3D Plot
    # if Plotvorticity == true
    #     omegaMag_obs, step_text_omega, fig_omega = Create_Plot3D(gridlengthX, gridlengthY, gridlengthZ, omegaMag; title="|ω|")
    #     screenOmega = GLMakie.Screen()
    #     display(screenOmega, fig_omega)
    # end

    if Plotvx==true
        #xy slice at z=midZ
        vx_xy_obs, step_text_vx_xy, fig_vx_xy = Create_Plot_XY(gridlengthX, gridlengthY, velocityX[:,:,midZ]; title="v_x at z=$(midZ)")
        screen_vx_xy = GLMakie.Screen()
        display(screen_vx_xy, fig_vx_xy)

        #xz slice at y=midY
        vx_xz_obs, step_text_vx_xz, fig_vx_xz = Create_Plot_XZ(gridlengthX, gridlengthZ, velocityX[:,midY,:]; title="v_x at y=$(midY)")
        screen_vx_xz = GLMakie.Screen()
        display(screen_vx_xz, fig_vx_xz)
    end

    if Plotvorticity == true
        #xy slice at z=midZ
        omega_xy_obs, step_text_omega_xy, fig_omega_xy = Create_Vorticity_XY(gridlengthX, gridlengthY, omegaMag[:,:,midZ]; title = "|ω| at z=$(midZ)")
        screen_omega_xy = GLMakie.Screen()
        display(screen_omega_xy, fig_omega_xy)

        #xz slice at y=midY
        omega_xz_obs, step_text_omega_xz, fig_omega_xz = Create_Vorticity_XZ(gridlengthX, gridlengthZ, omegaMag[:,midY,:]; title = "|ω| at y=$(midY)")
        screen_omega_xz = GLMakie.Screen()
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
    ###### NEW STABILIZATION ######
    for x in 2:gridlengthX-1 #Iteration über alle Zellen außer die Randzellen
        for y in 2:gridlengthY-1
            for z in 2:gridlengthZ-1

                # Compute macroscopic quantities
                rho[x,y,z] = f000[x,y,z] + 
                            (fm00[x,y,z] + fp00[x,y,z] + f0m0[x,y,z] + f0p0[x,y,z] + f00m[x,y,z] + f00p[x,y,z]) +
                            (fmm0[x,y,z] + fmp0[x,y,z] + fpm0[x,y,z] + fpp0[x,y,z] + 
                            fm0m[x,y,z] + fm0p[x,y,z] + fp0m[x,y,z] + fp0p[x,y,z] +
                            f0mm[x,y,z] + f0mp[x,y,z] + f0pm[x,y,z] + f0pp[x,y,z])
                
                u[x,y,z] = ((-fm00[x,y,z] + fp00[x,y,z]) +
                            (-fmm0[x,y,z] - fmp0[x,y,z] + fpm0[x,y,z] + fpp0[x,y,z]) +
                            (-fm0m[x,y,z] - fm0p[x,y,z] + fp0m[x,y,z] + fp0p[x,y,z])) / rho[x,y,z]
                
                v[x,y,z] = ((-f0m0[x,y,z] + f0p0[x,y,z]) +
                            (-fmm0[x,y,z] + fmp0[x,y,z] - fpm0[x,y,z] + fpp0[x,y,z]) +
                            (-f0mm[x,y,z] - f0mp[x,y,z] + f0pm[x,y,z] + f0pp[x,y,z])) / rho[x,y,z]
                
                w[x,y,z] = ((-f00m[x,y,z] + f00p[x,y,z]) +
                            (-fm0m[x,y,z] + fm0p[x,y,z] - fp0m[x,y,z] + fp0p[x,y,z]) +
                            (-f0mm[x,y,z] + f0mp[x,y,z] - f0pm[x,y,z] + f0pp[x,y,z])) / rho[x,y,z]
                
                # Pre-compute polynomial factors
                u2 = u[x,y,z] * u[x,y,z]
                v2 = v[x,y,z] * v[x,y,z]
                w2 = w[x,y,z] * w[x,y,z]
                
                Pm_u = 1 - 3*u[x,y,z] + 3*u2
                P0_u = -2 + 3*u2
                Pp_u = 1 + 3*u[x,y,z] + 3*u2
                
                Pm_v = 1 - 3*v[x,y,z] + 3*v2
                P0_v = -2 + 3*v2
                Pp_v = 1 + 3*v[x,y,z] + 3*v2
                
                Pm_w = 1 - 3*w[x,y,z] + 3*w2
                P0_w = -2 + 3*w2
                Pp_w = 1 + 3*w[x,y,z] + 3*w2
                
                # Push scheme: Rest particle (0,0,0) - weight 1/3
                f000S[x,y,z] = f000[x,y,z] + omega * (rho[x,y,z] * P0_u * P0_v * P0_w / 3.0 - f000[x,y,z])
                
                # Push scheme: Face neighbors - weight 1/18
                fm00S[x-1,y,z] = fm00[x,y,z] + omega * (rho[x,y,z] * Pm_u * P0_v * P0_w / 18.0 - fm00[x,y,z])
                fp00S[x+1,y,z] = fp00[x,y,z] + omega * (rho[x,y,z] * Pp_u * P0_v * P0_w / 18.0 - fp00[x,y,z])
                
                f0m0S[x,y-1,z] = f0m0[x,y,z] + omega * (rho[x,y,z] * P0_u * Pm_v * P0_w / 18.0 - f0m0[x,y,z])
                f0p0S[x,y+1,z] = f0p0[x,y,z] + omega * (rho[x,y,z] * P0_u * Pp_v * P0_w / 18.0 - f0p0[x,y,z])
                
                f00mS[x,y,z-1] = f00m[x,y,z] + omega * (rho[x,y,z] * P0_u * P0_v * Pm_w / 18.0 - f00m[x,y,z])
                f00pS[x,y,z+1] = f00p[x,y,z] + omega * (rho[x,y,z] * P0_u * P0_v * Pp_w / 18.0 - f00p[x,y,z])
                
                # Push scheme: Edge neighbors - weight 1/36
                # XY-plane edges
                fmm0S[x-1,y-1,z] = fmm0[x,y,z] + omega * (rho[x,y,z] * Pm_u * Pm_v * P0_w / 36.0 - fmm0[x,y,z])
                fmp0S[x-1,y+1,z] = fmp0[x,y,z] + omega * (rho[x,y,z] * Pm_u * Pp_v * P0_w / 36.0 - fmp0[x,y,z])
                fpm0S[x+1,y-1,z] = fpm0[x,y,z] + omega * (rho[x,y,z] * Pp_u * Pm_v * P0_w / 36.0 - fpm0[x,y,z])
                fpp0S[x+1,y+1,z] = fpp0[x,y,z] + omega * (rho[x,y,z] * Pp_u * Pp_v * P0_w / 36.0 - fpp0[x,y,z])
                
                # XZ-plane edges
                fm0mS[x-1,y,z-1] = fm0m[x,y,z] + omega * (rho[x,y,z] * Pm_u * P0_v * Pm_w / 36.0 - fm0m[x,y,z])
                fm0pS[x-1,y,z+1] = fm0p[x,y,z] + omega * (rho[x,y,z] * Pm_u * P0_v * Pp_w / 36.0 - fm0p[x,y,z])
                fp0mS[x+1,y,z-1] = fp0m[x,y,z] + omega * (rho[x,y,z] * Pp_u * P0_v * Pm_w / 36.0 - fp0m[x,y,z])
                fp0pS[x+1,y,z+1] = fp0p[x,y,z] + omega * (rho[x,y,z] * Pp_u * P0_v * Pp_w / 36.0 - fp0p[x,y,z])
                
                # YZ-plane edges
                f0mmS[x,y-1,z-1] = f0mm[x,y,z] + omega * (rho[x,y,z] * P0_u * Pm_v * Pm_w / 36.0 - f0mm[x,y,z])
                f0mpS[x,y-1,z+1] = f0mp[x,y,z] + omega * (rho[x,y,z] * P0_u * Pm_v * Pp_w / 36.0 - f0mp[x,y,z])
                f0pmS[x,y+1,z-1] = f0pm[x,y,z] + omega * (rho[x,y,z] * P0_u * Pp_v * Pm_w / 36.0 - f0pm[x,y,z])
                f0ppS[x,y+1,z+1] = f0pp[x,y,z] + omega * (rho[x,y,z] * P0_u * Pp_v * Pp_w / 36.0 - f0pp[x,y,z])
            end
        end
    end

    ###### (new) Boundary Conditions ######
    #Bounceback walls

    for x in 1:gridlengthX
        for y in 1:gridlengthY
            #bottom wall (z=1)
            f00pS[x,y,1] = f00mS[x,y,1] #mitte
            f0mpS[x,y,1] = f0pmS[x,y,1] #4 kanten oben -> unten
            fm0pS[x,y,1] = fp0mS[x,y,1] 
            f0ppS[x,y,1] = f0mmS[x,y,1]
            fp0pS[x,y,1] = fm0mS[x,y,1]
            
            #top wall (z=rows)
            f00mS[x,y,gridlengthZ] = f00pS[x,y,gridlengthZ]  #mitte
            f0pmS[x,y,gridlengthZ] = f0mpS[x,y,gridlengthZ]  #4 kanten unten -> oben
            fp0mS[x,y,gridlengthZ] = fm0pS[x,y,gridlengthZ]  
            f0mmS[x,y,gridlengthZ] = f0ppS[x,y,gridlengthZ]
            fm0mS[x,y,gridlengthZ] = fp0pS[x,y,gridlengthZ] 
            

        end
    end
    
    for x in 1:gridlengthX
        for z in 1:gridlengthZ
            #front wall (y=1)
            f0m0S[x,1,z] = f0p0S[x,1,z]   #mitte
            f0mpS[x,1,z] = f0pmS[x,1,z]   #4 kanten
            fpm0S[x,1,z] = fmp0S[x,1,z]
            f0mmS[x,1,z] = f0ppS[x,1,z]
            fmm0S[x,1,z] = fpp0S[x,1,z]

            #back wall (y=)
            f0p0S[x,gridlengthY,z] = f0m0S[x,gridlengthY,z]   #mitte
            f0pmS[x,gridlengthY,z] = f0mpS[x,gridlengthY,z]   #4 kanten
            fmp0S[x,gridlengthY,z] = fpm0S[x,gridlengthY,z]
            f0ppS[x,gridlengthY,z] = f0mmS[x,gridlengthY,z]
            fpp0S[x,gridlengthY,z] = fmm0S[x,gridlengthY,z]
        end
    end
    
    #Bounceback sphere
    for x in 1:gridlengthX
        for y in 1:gridlengthY
            for z in 1:gridlengthZ
                if sphere[x,y,z]
                    # Cell faces - swap
                    fp00S[x,y,z], fm00S[x,y,z] = fm00S[x,y,z], fp00S[x,y,z]  # X-Richtung
                    f0p0S[x,y,z], f0m0S[x,y,z] = f0m0S[x,y,z], f0p0S[x,y,z]  # Y-Richtung
                    f00pS[x,y,z], f00mS[x,y,z] = f00mS[x,y,z], f00pS[x,y,z]  # Z-Richtung
                    
                    # XY-plane edges
                    fpp0S[x,y,z], fmm0S[x,y,z] = fmm0S[x,y,z], fpp0S[x,y,z]   
                    fpm0S[x,y,z], fmp0S[x,y,z] = fmp0S[x,y,z], fpm0S[x,y,z]   
                    
                    # XZ-plane edges
                    fp0pS[x,y,z], fm0mS[x,y,z] = fm0mS[x,y,z], fp0pS[x,y,z]  
                    fp0mS[x,y,z], fm0pS[x,y,z] = fm0pS[x,y,z], fp0mS[x,y,z]  
                    
                    # YZ-plane edges
                    f0ppS[x,y,z], f0mmS[x,y,z] = f0mmS[x,y,z], f0ppS[x,y,z]  
                    f0pmS[x,y,z], f0mpS[x,y,z] = f0mpS[x,y,z], f0pmS[x,y,z] 
                end 

            end
        end
    end 
    ###### (new) Boundary Conditions ######

    #Swap copy f_eq to new distributions

    f000 .= f000S
    fm00 .= fm00S
    fp00 .= fp00S
    f0m0 .= f0m0S
    f0p0 .= f0p0S
    f00m .= f00mS
    f00p .= f00pS
    fmm0 .= fmm0S
    fmp0 .= fmp0S
    fpm0 .= fpm0S
    fpp0 .= fpp0S
    fm0m .= fm0mS
    fm0p .= fm0pS
    fp0m .= fp0mS
    fp0p .= fp0pS
    f0mm .= f0mmS
    f0mp .= f0mpS
    f0pm .= f0pmS
    f0pp .= f0ppS

    ###### NEW STABILIZATION #######

    # # Get Macroscopic values 3D
    # global densityGrid = sum(distributions, dims=4);
    # velocityX .= (1 ./ densityGrid) .* sum(distributions.*velocity_vector_x, dims=4); 
    # velocityY .= (1 ./ densityGrid) .* sum(distributions.*velocity_vector_y, dims=4); 
    # velocityZ .= (1 ./ densityGrid) .* sum(distributions.*velocity_vector_z, dims=4); 

    # #DEBUGGING############
    # if i % 20 == 0
    #     @views begin
    #         xy_min, xy_max = extrema(velocityX[:,:,midZ])
    #         xz_min, xz_max = extrema(velocityX[:,midY,:])
    #     end

    # println("vx XY min/max : ($(xy_min), $(xy_max)) | XZ min/max: ($(xz_min), $(xz_max)) at step $i")
    # end
    # #DEBUGGING############


    # ## Apply Collision 3D
    # # Compute equilibrium state
    # dotprod_velocities .= (velocity_vector_x .* velocityX) .+ (velocity_vector_y .* velocityY) .+ (velocity_vector_z .* velocityZ);
    # distributions_equilibrium .= weights .* densityGrid .*(1 .+ 3 .*dotprod_velocities .+ 4.5 .*dotprod_velocities.^2 .- 1.5 .*(velocityX.^2 .+ velocityY.^2 + velocityZ.^2));
    # # Relax towards equilibrium
    # distributions .+= -(1/τ) .* (distributions .- distributions_equilibrium);

    # # Stream 
    # for j in 1:Q
    #     distributions[:,:,:,j] = circshift(distributions[:,:,:,j], (velocity_vector_x[j], velocity_vector_y[j], velocity_vector_z[j]))
    # end

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

    # #No Slip Walls 3D
    # distributions[walls, 1:Q] .= distributions[walls, [1,3,2,5,4,7,6,9,8,11,10,13,12,15,14,17,16,19,18]];
    
    # # Apply object boundary condition 3D
    # distributions[sphere, 1:Q] .= distributions[sphere, [1,3,2,5,4,7,6,9,8,11,10,13,12,15,14,17,16,19,18]];

        # Plot of the field
    if ((i % 10 == 0)) || (i == simulationTime)

        #Copy velocities for plotting
        velocityX .= u
        velocityY .= v
        velocityZ .= w

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