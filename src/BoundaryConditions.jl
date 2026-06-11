module BoundaryConditions

export compute_object_boundary_data, apply_bouzidi_bc_3d!

function _ray_cylinder_q(lx0, ly0, dlx, dly, z0_phys, dz_phys, 
                         cylinder_radius, z_bot_phys, z_top_phys)
    
    q = Inf

    # ray: r(t) = l0 + t*dl
    # circle of cylinder: x^2 + y^2 = R^2
    # ray in circle equation:
    # (x0 + t*dx)^2 + (y0 + t*dy)^2 - R^2 = 0
    # a*t^2 + b*t + c = 0
    # check for intersections with discriminant

    a = dlx^2 + dly^2
    if a > 1e-10
        b = 2.0 * (lx0* dlx + ly0 * dly)
        c = lx0^2 + ly0^2 - cylinder_radius^2

        disc = b^2 - 4.0 * a * c

        if disc >= 0.0
            sqrt_disc = sqrt(disc)
            t1 = (-b - sqrt_disc) / (2 * a)
            t2 = (-b + sqrt_disc) / (2 * a)

            for t in (t1, t2)
                if 0.0 < t <= 1.0 + 1e-10
                    z_at_t = z0_phys + t * dz_phys
                    if z_at_t >= z_bot_phys - 1e-10 && z_at_t <= z_top_phys + 1e-10
                        q = min(q, t)
                    end
                end
            end
        end
    end
    
    # ray: z(t) = z0 + t * dz
    # z_cap: z-coord for bot/top
    # z_cap = z0 + t * dz <-> t = (z_cap - z0_phys) / dz_phys
    # then: check if x/y are INSIDE of circle area
    if abs(dz_phys) > 1e-10
        for z_cap in (z_bot_phys, z_top_phys)
            t = (z_cap - z0_phys) / dz_phys
            if 0.0 < t <= 1.0 + 1e-10
               x_at_t = lx0 + t * dlx
               y_at_t = ly0 + t * dly
               if x_at_t^2 + y_at_t^2 <= cylinder_radius^2 + 1e-10
                    q  = min(q, t)
               end
            end
        end
    end

    return q    
end

function compute_object_boundary_data(gridlengthX, gridlengthY, gridlengthZ,
                                      cylinder_x, cylinder_y, cylinder_radius,
                                      cylinder_z_bot, cylinder_z_top, is_object, delta_x)


    # boundary_data = []
    boundary_data = Tuple{Int,Int,Int,Int,Int,Int,Int,Int,Float64}[]

    # D3Q19 - without (0,0,0)
    directions = (
    ( 1, 0, 0,  2,  1), (-1, 0, 0,  1,  2),
    ( 0, 1, 0,  4,  3), ( 0,-1, 0,  3,  4),
    ( 0, 0, 1,  6,  5), ( 0, 0,-1,  5,  6),
    # Edges
    (-1,-1, 0,  7, 10), (-1, 1, 0,  8,  9),
    ( 1,-1, 0,  9,  8), ( 1, 1, 0, 10,  7),
    (-1, 0,-1, 11, 14), (-1, 0, 1, 12, 13),
    ( 1, 0,-1, 13, 12), ( 1, 0, 1, 14, 11),
    ( 0,-1,-1, 15, 18), ( 0,-1, 1, 16, 17),
    ( 0, 1,-1, 17, 16), ( 0, 1, 1, 18, 15)
    )
    
    # Loop over fluid nodes
    @inbounds for x in 2:gridlengthX-1, y in 2:gridlengthY-1, z in 2:gridlengthZ-1

        # Skip if in cylinder
        if is_object[x, y, z]
            continue
        end
   
        # convert to physical coordinates
        x0_phys = (x-2) * delta_x
        y0_phys = (y-2) * delta_x
        z0_phys = (z-2) * delta_x

        # push so that cylinder middle = (0,0)

        lx0 = x0_phys - cylinder_x
        ly0 = y0_phys - cylinder_y

        for (cx, cy, cz, idx_toward, idx_reflect) in directions

            # neighbor positions
            nx = x + cx
            ny = y + cy
            nz = z + cz

            # Bounds check: is neighbor inside grid
            if nx < 1 || nx > gridlengthX || ny < 1 || ny > gridlengthY || nz < 1 || nz > gridlengthZ
                continue
            end

            # Is neighbor in cylinder?
            if !is_object[nx, ny, nz]
                continue
            end

            # neighbor position physical
            x1_phys = (nx-2) * delta_x
            y1_phys = (ny-2) * delta_x
            z1_phys = (nz-2) * delta_x

            # local coordinates of neighbor
            lx1 = x1_phys - cylinder_x
            ly1 = y1_phys - cylinder_y
            dz_phys = z1_phys - z0_phys

            # Ray-casting -> q
            q = _ray_cylinder_q(lx0, ly0, lx1-lx0, ly1-ly0, 
                                z0_phys, dz_phys, 
                                cylinder_radius, cylinder_z_bot, cylinder_z_top)


            # for top/bot of cylinder and edge:
            # ray doenst cross curved surface of cylinder, but we are at fluid node and neigbor in direction is solid!
            # so "normal" bb -> q = 0.5
            if q == Inf
                q = 0.5
            end
            
            # check valid q then push to boundary_data (write to end of array)
            if 0.0 < q <= 1.0+1e-10
                push!(boundary_data, (x, y, z, idx_toward, idx_reflect, cx, cy, cz, q))
            end
        end #for (cx, cy, cz)
    end #@inbounds for x,y,z

    println("✓ Bouzidi BC: $(length(boundary_data)) boundary nodes found")
    return boundary_data
end


# FAST VERSION?
function apply_bouzidi_bc_3d!(boundary_data, fS)

    F_x = 0.0
    F_y = 0.0

      @inbounds for (x, y, z, idx_toward, idx_reflect, cx, cy, cz, q) in boundary_data
        q_toward_solid = idx_toward  + 1   # 1-based non-rest index -> Q index in fS
        q_reflected    = idx_reflect + 1

        q2 = 2.0 * q
        f_at_solid = fS[q_toward_solid, x + cx, y + cy, z + cz]

        if q < 0.5
            f_new = q2 * f_at_solid + (1.0 - q2) * fS[q_toward_solid, x, y, z]
        else
            iq2 = 1.0 / q2
            r = (q2 - 1.0) * iq2
            f_opposite = fS[q_reflected, x - cx, y - cy, z - cz]
            f_new = iq2 * f_at_solid + r * f_opposite
        end

        fS[q_reflected, x, y, z] = f_new

        # # Force calculation
        # F_x += cx * (fS[q_toward_solid, x, y, z] + f_new)
        # F_y += cy * (fS[q_toward_solid, x, y, z] + f_new)

        F_x += cx * (f_at_solid + f_new)
        F_y += cy * (f_at_solid + f_new)



    end

    return F_x,  F_y
end




# SLOW VERSION?
# function apply_bouzidi_bc_3d!(boundary_data, f000S, fm00S, fp00S, f0m0S, f0p0S, f00mS, f00pS,
#                               fmm0S, fmp0S, fpm0S, fpp0S,
#                               fm0mS, fm0pS, fp0mS, fp0pS,
#                               f0mmS, f0mpS, f0pmS, f0ppS)

#     @inbounds for (x, y, z, cx, cy, cz, q) in boundary_data
        
#         q2 = 2.0 * q2
#         f_in = 0.0
#         f_out = 0.0
    
#         if q < 0.5
#             # q < 0.5: f_out = 2q*f_in + (1-2q)*f_fluid
            
#             if cx == 1 && cy == 0 && cz == 0
#                 f_in = fp00S[x+1, y, z]
#                 f_out = q2 * f_in + (1.0 - q2) * fp00S[x, y, z]
#                 fm00S[x, y, z] = f_out
                
#             elseif cx == -1 && cy == 0 && cz == 0
#                 f_in = fm00S[x-1, y, z]
#                 f_out = q2 * f_in + (1.0 - q2) * fm00S[x, y, z]
#                 fp00S[x, y, z] = f_out
                
#             elseif cx == 0 && cy == 1 && cz == 0
#                 f_in = f0p0S[x, y+1, z]
#                 f_out = q2 * f_in + (1.0 - q2) * f0p0S[x, y, z]
#                 f0m0S[x, y, z] = f_out
                
#             elseif cx == 0 && cy == -1 && cz == 0
#                 f_in = f0m0S[x, y-1, z]
#                 f_out = q2 * f_in + (1.0 - q2) * f0m0S[x, y, z]
#                 f0p0S[x, y, z] = f_out
                
#             elseif cx == 0 && cy == 0 && cz == 1
#                 f_in = f00pS[x, y, z+1]
#                 f_out = q2 * f_in + (1.0 - q2) * f00pS[x, y, z]
#                 f00mS[x, y, z] = f_out
                
#             elseif cx == 0 && cy == 0 && cz == -1
#                 f_in = f00mS[x, y, z-1]
#                 f_out = q2 * f_in + (1.0 - q2) * f00mS[x, y, z]
#                 f00pS[x, y, z] = f_out
                
#             # XY edges
#             elseif cx == 1 && cy == 1 && cz == 0
#                 f_in = fpp0S[x+1, y+1, z]
#                 f_out = q2 * f_in + (1.0 - q2) * fpp0S[x, y, z]
#                 fmm0S[x, y, z] = f_out
                
#             elseif cx == -1 && cy == -1 && cz == 0
#                 f_in = fmm0S[x-1, y-1, z]
#                 f_out = q2 * f_in + (1.0 - q2) * fmm0S[x, y, z]
#                 fpp0S[x, y, z] = f_out
                
#             elseif cx == 1 && cy == -1 && cz == 0
#                 f_in = fpm0S[x+1, y-1, z]
#                 f_out = q2 * f_in + (1.0 - q2) * fpm0S[x, y, z]
#                 fmp0S[x, y, z] = f_out
                
#             elseif cx == -1 && cy == 1 && cz == 0
#                 f_in = fmp0S[x-1, y+1, z]
#                 f_out = q2 * f_in + (1.0 - q2) * fmp0S[x, y, z]
#                 fpm0S[x, y, z] = f_out
                
#             # XZ edges
#             elseif cx == 1 && cy == 0 && cz == 1
#                 f_in = fp0pS[x+1, y, z+1]
#                 f_out = q2 * f_in + (1.0 - q2) * fp0pS[x, y, z]
#                 fm0mS[x, y, z] = f_out
                
#             elseif cx == -1 && cy == 0 && cz == -1
#                 f_in = fm0mS[x-1, y, z-1]
#                 f_out = q2 * f_in + (1.0 - q2) * fm0mS[x, y, z]
#                 fp0pS[x, y, z] = f_out
                
#             elseif cx == 1 && cy == 0 && cz == -1
#                 f_in = fp0mS[x+1, y, z-1]
#                 f_out = q2 * f_in + (1.0 - q2) * fp0mS[x, y, z]
#                 fm0pS[x, y, z] = f_out
                
#             elseif cx == -1 && cy == 0 && cz == 1
#                 f_in = fm0pS[x-1, y, z+1]
#                 f_out = q2 * f_in + (1.0 - q2) * fm0pS[x, y, z]
#                 fp0mS[x, y, z] = f_out
                
#             # YZ edges
#             elseif cx == 0 && cy == 1 && cz == 1
#                 f_in = f0ppS[x, y+1, z+1]
#                 f_out = q2 * f_in + (1.0 - q2) * f0ppS[x, y, z]
#                 f0mmS[x, y, z] = f_out
                
#             elseif cx == 0 && cy == -1 && cz == -1
#                 f_in = f0mmS[x, y-1, z-1]
#                 f_out = q2 * f_in + (1.0 - q2) * f0mmS[x, y, z]
#                 f0ppS[x, y, z] = f_out
                
#             elseif cx == 0 && cy == 1 && cz == -1
#                 f_in = f0pmS[x, y+1, z-1]
#                 f_out = q2 * f_in + (1.0 - q2) * f0pmS[x, y, z]
#                 f0mpS[x, y, z] = f_out
                
#             elseif cx == 0 && cy == -1 && cz == 1
#                 f_in = f0mpS[x, y-1, z+1]
#                 f_out = q2 * f_in + (1.0 - q2) * f0mpS[x, y, z]
#                 f0pmS[x, y, z] = f_out
#             end
            
#         else
#             # q >= 0.5: f_out = f_in/(2q) + (2q-1)/(2q)*f_opposite
#             iq2 = 1.0 / q2
#             r = (q2 - 1.0) * iq2
            
#             if cx == 1 && cy == 0 && cz == 0
#                 f_in = fp00S[x+1, y, z]
#                 f_out = iq2 * f_in + r * fm00S[x-1, y, z]
#                 fm00S[x, y, z] = f_out
                
#             elseif cx == -1 && cy == 0 && cz == 0
#                 f_in = fm00S[x-1, y, z]
#                 f_out = iq2 * f_in + r * fp00S[x+1, y, z]
#                 fp00S[x, y, z] = f_out
                
#             elseif cx == 0 && cy == 1 && cz == 0
#                 f_in = f0p0S[x, y+1, z]
#                 f_out = iq2 * f_in + r * f0m0S[x, y-1, z]
#                 f0m0S[x, y, z] = f_out
                
#             elseif cx == 0 && cy == -1 && cz == 0
#                 f_in = f0m0S[x, y-1, z]
#                 f_out = iq2 * f_in + r * f0p0S[x, y+1, z]
#                 f0p0S[x, y, z] = f_out
                
#             elseif cx == 0 && cy == 0 && cz == 1
#                 f_in = f00pS[x, y, z+1]
#                 f_out = iq2 * f_in + r * f00mS[x, y, z-1]
#                 f00mS[x, y, z] = f_out
                
#             elseif cx == 0 && cy == 0 && cz == -1
#                 f_in = f00mS[x, y, z-1]
#                 f_out = iq2 * f_in + r * f00pS[x, y, z+1]
#                 f00pS[x, y, z] = f_out
                
#             # XY edges
#             elseif cx == 1 && cy == 1 && cz == 0
#                 f_in = fpp0S[x+1, y+1, z]
#                 f_out = iq2 * f_in + r * fmm0S[x-1, y-1, z]
#                 fmm0S[x, y, z] = f_out
                
#             elseif cx == -1 && cy == -1 && cz == 0
#                 f_in = fmm0S[x-1, y-1, z]
#                 f_out = iq2 * f_in + r * fpp0S[x+1, y+1, z]
#                 fpp0S[x, y, z] = f_out
                
#             elseif cx == 1 && cy == -1 && cz == 0
#                 f_in = fpm0S[x+1, y-1, z]
#                 f_out = iq2 * f_in + r * fmp0S[x-1, y+1, z]
#                 fmp0S[x, y, z] = f_out
                
#             elseif cx == -1 && cy == 1 && cz == 0
#                 f_in = fmp0S[x-1, y+1, z]
#                 f_out = iq2 * f_in + r * fpm0S[x+1, y-1, z]
#                 fpm0S[x, y, z] = f_out
                
#             # XZ edges
#             elseif cx == 1 && cy == 0 && cz == 1
#                 f_in = fp0pS[x+1, y, z+1]
#                 f_out = iq2 * f_in + r * fm0mS[x-1, y, z-1]
#                 fm0mS[x, y, z] = f_out
                
#             elseif cx == -1 && cy == 0 && cz == -1
#                 f_in = fm0mS[x-1, y, z-1]
#                 f_out = iq2 * f_in + r * fp0pS[x+1, y, z+1]
#                 fp0pS[x, y, z] = f_out
                
#             elseif cx == 1 && cy == 0 && cz == -1
#                 f_in = fp0mS[x+1, y, z-1]
#                 f_out = iq2 * f_in + r * fm0pS[x-1, y, z+1]
#                 fm0pS[x, y, z] = f_out
                
#             elseif cx == -1 && cy == 0 && cz == 1
#                 f_in = fm0pS[x-1, y, z+1]
#                 f_out = iq2 * f_in + r * fp0mS[x+1, y, z-1]
#                 fp0mS[x, y, z] = f_out
                
#             # YZ edges
#             elseif cx == 0 && cy == 1 && cz == 1
#                 f_in = f0ppS[x, y+1, z+1]
#                 f_out = iq2 * f_in + r * f0mmS[x, y-1, z-1]
#                 f0mmS[x, y, z] = f_out
                
#             elseif cx == 0 && cy == -1 && cz == -1
#                 f_in = f0mmS[x, y-1, z-1]
#                 f_out = iq2 * f_in + r * f0ppS[x, y+1, z+1]
#                 f0ppS[x, y, z] = f_out
                
#             elseif cx == 0 && cy == 1 && cz == -1
#                 f_in = f0pmS[x, y+1, z-1]
#                 f_out = iq2 * f_in + r * f0mpS[x, y-1, z+1]
#                 f0mpS[x, y, z] = f_out
                
#             elseif cx == 0 && cy == -1 && cz == 1
#                 f_in = f0mpS[x, y-1, z+1]
#                 f_out = iq2 * f_in + r * f0pmS[x, y+1, z-1]
#                 f0pmS[x, y, z] = f_out
#             end
#         end
#     end
    
#     return nothing

end #module BoundaryConditions