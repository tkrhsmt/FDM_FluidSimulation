# --------------------------------------------------
# MODULE    Log
# julia code
#
# DATE : 2024/8/23
# --------------------------------------------------

module Log

export init_log, calculation_log, output_log, final_log

# 2nd-order central difference
ddx1_1(u1, u2, dx) = (u2 - u1) / (dx)

function init_log(prm)

    println("==================================================")
    println("                 FLOW SIMULATION")
    println("==================================================")

    println("SYSTEM   :: CPU")
    println("EQUATION :: Navier-Stokes Equation")
    println("SOLVER   :: 2nd order central difference")

    println("--------------------------------------------------")

    println("                   PARAMETER")
    println("nx  :: $(prm.nx)")
    println("ny  :: $(prm.ny)")
    println("nz  :: $(prm.nz)\n")
    println("lx  :: $(prm.lx)")
    println("ly  :: $(prm.ly)")
    println("lz  :: $(prm.lz)\n")

    println("dt  :: $(prm.dt)")
    println("rep :: $(prm.istart) -- $(prm.iend)")
    println("re  :: $(prm.re)")

    println("--------------------------------------------------")
end

function calculation_log(ux, uy, uz, pp, prm, rep, time)

    println("$time / $(prm.iend)")

    # spatial step size
    dx = prm.lx / prm.nx
    dy = prm.ly / prm.ny
    dz = prm.lz / prm.nz

    # the range of pressure components not including the boundaries
    nx1 = 2
    nx2 = prm.nx + 1
    ny1 = 2
    ny2 = prm.ny + 1
    nz1 = 2
    nz2 = prm.nz + 1

    # maximum velocity divergence
    div_max = 0.0
    # average velocity divergence
    div_ave = 0.0

    for k in nz1:nz2
        for j in ny1:ny2
            for i in nx1:nx2
                # dudx
                tmp1 = ddx1_1(ux[i, j, k], ux[i+1, j, k], dx)
                # dudx + dvdy
                tmp1 += ddx1_1(uy[i, j, k], uy[i, j+1, k], dy)
                # dudx + dvdy + dwdz
                tmp1 += ddx1_1(uz[i, j, k], uz[i, j, k+1], dz)

                if div_max < tmp1
                    div_max = tmp1
                end
                div_ave += tmp1
            end
        end
    end
    div_ave /= (prm.nx * prm.ny * prm.nz)

    println("div mean   :: $div_ave")
    println("div max    :: $div_max")

    cfl_x = maximum(ux[nx1:nx2, ny1:ny2, nz1:nz2] / dx * prm.dt)
    println("CFL x      :: $cfl_x")
    cfl_y = maximum(uy[nx1:nx2, ny1:ny2, nz1:nz2] / dy * prm.dt)
    println("CFL y      :: $cfl_y")
    cfl_z = maximum(uz[nx1:nx2, ny1:ny2, nz1:nz2] / dz * prm.dt)
    println("CFL z      :: $cfl_z")

    @assert cfl_x < 1.0 "ERROR : CFL X IS OVERED!"
    @assert cfl_y < 1.0 "ERROR : CFL Y IS OVERED!"
    @assert cfl_z < 1.0 "ERROR : CFL Z IS OVERED!"

    ux_max = maximum(ux[nx1:nx2, ny1:ny2, nz1:nz2])
    uy_max = maximum(uy[nx1:nx2, ny1:ny2, nz1:nz2])
    uz_max = maximum(uz[nx1:nx2, ny1:ny2, nz1:nz2])
    pp_max = maximum(pp[nx1:nx2, ny1:ny2, nz1:nz2])
    ux_min = minimum(ux[nx1:nx2, ny1:ny2, nz1:nz2])
    uy_min = minimum(uy[nx1:nx2, ny1:ny2, nz1:nz2])
    uz_min = minimum(uz[nx1:nx2, ny1:ny2, nz1:nz2])
    pp_min = minimum(pp[nx1:nx2, ny1:ny2, nz1:nz2])

    println("ux min max :: $ux_min\t$ux_max")
    println("uy min max :: $uy_min\t$uy_max")
    println("uz min max :: $uz_min\t$uz_max")
    println("pp min max :: $pp_min\t$pp_max")

    println("repeat num :: $rep")


    println("--------------------------------------------------")

end

function output_log(time, prm)

    println("FILE OUTPUT SUCCESSED!")
    println("$time -> /data/output-$(Int(time / prm.ioutput))")
    println("--------------------------------------------------")

end

function final_log()

    println("                   PROGRAM END")
    println("==================================================")

end


end
