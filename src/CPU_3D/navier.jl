# --------------------------------------------------
# MODULE    Navier
# julia code
#
# DATE : 2024/8/23
# --------------------------------------------------

module Navier

export first_velocity!, second_velocity!

include("les_model.jl")
using .Les

# 2nd-order central difference
ddx1_1(u1, u2, dx) = (u2 - u1) / (dx)
ddx2_1(u1, u2, u3, dx) = (u3 - 2 * u2 + u1) / (dx^2)

# 2nd-order central interpolation
int1(u1, u2) = (u1 + u2) / 2.0

function first_velocity!(ux1, uy1, uz1, ux2, uy2, uz2, prm, boundary_ux, boundary_uy, boundary_uz, force_ux, force_uy, force_uz)

    # the range of x-direction velocity components not including the boundaries
    nx1 = 2
    nx2 = prm.nx + 2
    ny1 = 2
    ny2 = prm.ny + 1
    nz1 = 2
    nz2 = prm.nz + 1

    # spatial step size
    dx = prm.lx / prm.nx
    dy = prm.ly / prm.ny
    dz = prm.lz / prm.nz

    # x-direction pressureless velocity
    for k in nz1:nz2
        for j in ny1:ny2
            for i in nx1:nx2

                # first advection term
                tmp1 = int1(ux1[i, j, k], ux1[i+1, j, k]) * ddx1_1(ux1[i, j, k], ux1[i+1, j, k], dx)
                tmp2 = int1(ux1[i-1, j, k], ux1[i, j, k]) * ddx1_1(ux1[i-1, j, k], ux1[i, j, k], dx)
                ADVX = int1(tmp1, tmp2)

                # second advection term
                tmp1 = int1(uy1[i-1, j, k], uy1[i, j, k]) * ddx1_1(ux1[i, j-1, k], ux1[i, j, k], dy)
                tmp2 = int1(uy1[i-1, j+1, k], uy1[i, j+1, k]) * ddx1_1(ux1[i, j, k], ux1[i, j+1, k], dy)
                ADVY = int1(tmp1, tmp2)

                # third advection term
                tmp1 = int1(uz1[i-1, j, k], uz1[i, j, k]) * ddx1_1(ux1[i, j, k-1], ux1[i, j, k], dz)
                tmp2 = int1(uz1[i-1, j, k+1], uz1[i, j, k+1]) * ddx1_1(ux1[i, j, k], ux1[i, j, k+1], dz)
                ADVZ = int1(tmp1, tmp2)

                # first diffusion term
                DIFX = ddx2_1(ux1[i+1, j, k], ux1[i, j, k], ux1[i-1, j, k], dx)

                # second diffusion term
                DIFY = ddx2_1(ux1[i, j+1, k], ux1[i, j, k], ux1[i, j-1, k], dy)

                # third diffusion term
                DIFZ = ddx2_1(ux1[i, j, k+1], ux1[i, j, k], ux1[i, j, k-1], dz)

                # next step pressureless x-direction velocity
                ux2[i, j, k] = ux1[i, j, k] - prm.dt * (ADVX + ADVY + ADVZ - 1 / prm.re * (DIFX + DIFY + DIFZ))

            end
        end
    end

    # add x-direction velocity force
    force_ux(ux2, nx1, nx2, ny1, ny2, nz1, nz2, prm.dt)

    # the range of y-direction velocity components not including the boundaries
    nx1 = 2
    nx2 = prm.nx + 1
    ny1 = 2
    ny2 = prm.ny + 2
    nz1 = 2
    nz2 = prm.nz + 1

    # y-direction pressureless velocity
    for k in nz1:nz2
        for j in ny1:ny2
            for i in nx1:nx2

                # first advection term
                tmp1 = int1(ux1[i, j-1, k], ux1[i, j, k]) * ddx1_1(uy1[i-1, j, k], uy1[i, j, k], dx)
                tmp2 = int1(ux1[i+1, j-1, k], ux1[i+1, j, k]) * ddx1_1(uy1[i, j, k], uy1[i+1, j, k], dx)
                ADVX = int1(tmp1, tmp2)

                # second advection term
                tmp1 = int1(uy1[i, j-1, k], uy1[i, j, k]) * ddx1_1(uy1[i, j-1, k], uy1[i, j, k], dy)
                tmp2 = int1(uy1[i, j, k], uy1[i, j+1, k]) * ddx1_1(uy1[i, j, k], uy1[i, j+1, k], dy)
                ADVY = int1(tmp1, tmp2)

                # third advection term
                tmp1 = int1(uz1[i, j-1, k], uz1[i, j, k]) * ddx1_1(uy1[i, j, k-1], uy1[i, j, k], dz)
                tmp2 = int1(uz1[i, j-1, k+1], uz1[i, j, k+1]) * ddx1_1(uy1[i, j, k], uy1[i, j, k+1], dz)
                ADVZ = int1(tmp1, tmp2)

                # first diffusion term
                DIFX = ddx2_1(uy1[i+1, j, k], uy1[i, j, k], uy1[i-1, j, k], dx)

                # second diffusion term
                DIFY = ddx2_1(uy1[i, j+1, k], uy1[i, j, k], uy1[i, j-1, k], dy)

                # third diffusion term
                DIFZ = ddx2_1(uy1[i, j, k+1], uy1[i, j, k], uy1[i, j, k-1], dx)

                # next step pressureless x-direction velocity
                uy2[i, j, k] = uy1[i, j, k] - prm.dt * (ADVX + ADVY + ADVZ - 1 / prm.re * (DIFX + DIFY + DIFZ))

            end
        end
    end

    # add y-direction velocity force
    force_uy(uy2, nx1, nx2, ny1, ny2, nz1, nz2, prm.dt)

    # the range of z-direction velocity components not including the boundaries
    nx1 = 2
    nx2 = prm.nx + 1
    ny1 = 2
    ny2 = prm.ny + 1
    nz1 = 2
    nz2 = prm.nz + 2


    # z-direction pressureless velocity
    for k in nz1:nz2
        for j in ny1:ny2
            for i in nx1:nx2

                # first advection term
                tmp1 = int1(ux1[i, j, k-1], ux1[i, j, k]) * ddx1_1(uz1[i-1, j, k], uz1[i, j, k], dx)
                tmp2 = int1(ux1[i+1, j, k-1], ux1[i+1, j, k]) * ddx1_1(uz1[i, j, k], uz1[i+1, j, k], dx)
                ADVX = int1(tmp1, tmp2)

                # second advection term
                tmp1 = int1(uy1[i, j, k-1], uy1[i, j, k]) * ddx1_1(uz1[i, j-1, k], uz1[i, j, k], dy)
                tmp2 = int1(uy1[i, j+1, k-1], uy1[i, j+1, k]) * ddx1_1(uz1[i, j, k], uz1[i, j+1, k], dy)
                ADVY = int1(tmp1, tmp2)

                # third advection term
                tmp1 = int1(uz1[i, j, k], uz1[i, j, k+1]) * ddx1_1(uz1[i, j, k], uz1[i, j, k+1], dz)
                tmp2 = int1(uz1[i, j, k-1], uz1[i, j, k]) * ddx1_1(uz1[i, j, k-1], uz1[i, j, k], dz)
                ADVZ = int1(tmp1, tmp2)

                # first diffusion term
                DIFX = ddx2_1(uz1[i+1, j, k], uz1[i, j, k], uz1[i-1, j, k], dx)

                # second diffusion term
                DIFY = ddx2_1(uz1[i, j+1, k], uz1[i, j, k], uz1[i, j-1, k], dy)

                # third diffusion term
                DIFZ = ddx2_1(uz1[i, j, k+1], uz1[i, j, k], uz1[i, j, k-1], dz)

                # next step pressureless x-direction velocity
                uz2[i, j, k] = uz1[i, j, k] - prm.dt * (ADVX + ADVY + ADVZ - 1 / prm.re * (DIFX + DIFY + DIFZ))

            end
        end
    end

    # add z-direction velocity force
    force_uz(uz2, nx1, nx2, ny1, ny2, nz1, nz2, prm.dt)

    # les model
    if prm.les == true
        dτ1, dτ2, dτ3 = les_model(ux1, uy1, uz1, prm)
        ux2 = ux2 - prm.dt * dτ1
        uy2 = uy2 - prm.dt * dτ2
        uz2 = uz2 - prm.dt * dτ3
    end

    # add x-direction boundary condition
    boundary_ux(ux2, nx1, nx2, ny1, ny2, nz1, nz2)

    # add y-direction boundary condition
    boundary_uy(uy2, nx1, nx2, ny1, ny2, nz1, nz2)

    # add z-direction boundary condition
    boundary_uz(uz2, nx1, nx2, ny1, ny2, nz1, nz2)

end

function second_velocity!(ux1, uy1, uz1, ux2, uy2, uz2, pp, prm, boundary_ux, boundary_uy, boundary_uz)

    # the range of x-direction velocity components not including the boundaries
    nx1 = 2
    nx2 = prm.nx + 2
    ny1 = 2
    ny2 = prm.ny + 1
    nz1 = 2
    nz2 = prm.nz + 1

    # spatial step size
    dx = prm.lx / prm.nx
    dy = prm.ly / prm.ny
    dz = prm.lz / prm.nz

    # x-direction pressureless velocity
    for k in nz1:nz2
        for j in ny1:ny2
            for i in nx1:nx2
                # calculate dpdx
                dpdx = ddx1_1(pp[i-1, j, k], pp[i, j, k], dx)
                ux2[i, j, k] = ux1[i, j, k] - prm.dt * dpdx
            end
        end
    end

    # add x-direction boundary condition
    boundary_ux(ux2, nx1, nx2, ny1, ny2, nz1, nz2)

    # the range of y-direction velocity components not including the boundaries
    nx1 = 2
    nx2 = prm.nx + 1
    ny1 = 2
    ny2 = prm.ny + 2
    nz1 = 2
    nz2 = prm.nz + 1

    # y-direction pressureless velocity
    for k in nz1:nz2
        for j in ny1:ny2
            for i in nx1:nx2
                # calculate dpdx
                dpdy = ddx1_1(pp[i, j-1, k], pp[i, j, k], dy)
                uy2[i, j, k] = uy1[i, j, k] - prm.dt * dpdy
            end
        end
    end

    # add y-direction boundary condition
    boundary_uy(uy2, nx1, nx2, ny1, ny2, nz1, nz2)

    # the range of z-direction velocity components not including the boundaries
    nx1 = 2
    nx2 = prm.nx + 1
    ny1 = 2
    ny2 = prm.ny + 1
    nz1 = 2
    nz2 = prm.nz + 2

    # z-direction pressureless velocity
    for k in nz1:nz2
        for j in ny1:ny2
            for i in nx1:nx2
                # calculate dpdx
                dpdz = ddx1_1(pp[i, j, k-1], pp[i, j, k], dz)
                uz2[i, j, k] = uz1[i, j, k] - prm.dt * dpdz
            end
        end
    end

    # add z-direction boundary condition
    boundary_uz(uz2, nx1, nx2, ny1, ny2, nz1, nz2)

end


end
