# --------------------------------------------------
# MODULE    Navier
# julia code
#
# DATE : 2024/8/23
# --------------------------------------------------

module Navier

include("input_param.jl")
using .Param

export first_velocity!, second_velocity!

# 2nd-order central difference
ddx1_1(u1, u2, dx) = (u2 - u1) / (dx)
ddx2_1(u1, u2, u3, dx) = (u3 - 2 * u2 + u1) / (dx^2)

# 2nd-order central interpolation
int1(u1, u2) = (u1 + u2) / 2.0

function first_velocity!(ux1, uy1, ux2, uy2, prm)

    # the range of x-direction velocity components not including the boundaries
    nx1 = 2
    nx2 = prm.nx + 2
    ny1 = 2
    ny2 = prm.ny + 1

    # spatial step size
    dx = prm.lx / prm.nx
    dy = prm.ly / prm.ny

    # x-direction pressureless velocity
    for j in ny1:ny2
        for i in nx1:nx2

            # first advection term
            tmp1 = int1(ux1[i, j], ux1[i+1, j]) * ddx1_1(ux1[i, j], ux1[i+1, j], dx)
            tmp2 = int1(ux1[i-1, j], ux1[i, j]) * ddx1_1(ux1[i-1, j], ux1[i, j], dx)
            ADVX = int1(tmp1, tmp2)

            # second advection term
            tmp1 = int1(uy1[i-1, j], uy1[i, j]) * ddx1_1(ux1[i, j-1], ux1[i, j], dy)
            tmp2 = int1(uy1[i-1, j+1], uy1[i, j+1]) * ddx1_1(ux1[i, j], ux1[i, j+1], dy)
            ADVY = int1(tmp1, tmp2)

            # first diffusion term
            DIFX = ddx2_1(ux1[i+1, j], ux1[i, j], ux1[i-1, j], dx)

            # second diffusion term
            DIFY = ddx2_1(ux1[i, j+1], ux1[i, j], ux1[i, j-1], dy)

            # next step pressureless x-direction velocity
            ux2[i, j] = ux1[i, j] - prm.dt * (ADVX + ADVY - 1 / prm.re * (DIFX + DIFY))

        end
    end

    # add x-direction boundary condition
    boundary_ux(ux2, nx1, nx2, ny1, ny2)

    # the range of y-direction velocity components not including the boundaries
    nx1 = 2
    nx2 = prm.nx + 1
    ny1 = 2
    ny2 = prm.ny + 2

    # y-direction pressureless velocity
    for j in ny1:ny2
        for i in nx1:nx2

            # first advection term
            tmp1 = int1(ux1[i, j-1], ux1[i, j]) * ddx1_1(uy1[i-1, j], uy1[i, j], dx)
            tmp2 = int1(ux1[i+1, j-1], ux1[i+1, j]) * ddx1_1(uy1[i, j], uy1[i+1, j], dx)
            ADVX = int1(tmp1, tmp2)

            # second advection term
            tmp1 = int1(uy1[i, j-1], uy1[i, j]) * ddx1_1(uy1[i, j-1], uy1[i, j], dy)
            tmp2 = int1(uy1[i, j], uy1[i, j+1]) * ddx1_1(uy1[i, j], uy1[i, j+1], dy)
            ADVY = int1(tmp1, tmp2)

            # first diffusion term
            DIFX = ddx2_1(uy1[i+1, j], uy1[i, j], uy1[i-1, j], dx)

            # second diffusion term
            DIFY = ddx2_1(uy1[i, j+1], uy1[i, j], uy1[i, j-1], dy)

            # next step pressureless x-direction velocity
            uy2[i, j] = uy1[i, j] - prm.dt * (ADVX + ADVY - 1 / prm.re * (DIFX + DIFY))

        end
    end

    # add y-direction boundary condition
    boundary_uy(uy2, nx1, nx2, ny1, ny2)

end

function second_velocity!(ux1, uy1, ux2, uy2, pp, prm)

    # the range of x-direction velocity components not including the boundaries
    nx1 = 2
    nx2 = prm.nx + 2
    ny1 = 2
    ny2 = prm.ny + 1

    # spatial step size
    dx = prm.lx / prm.nx
    dy = prm.ly / prm.ny

    # x-direction pressureless velocity
    for j in ny1:ny2
        for i in nx1:nx2
            # calculate dpdx
            dpdx = ddx1_1(pp[i-1, j], pp[i, j], dx)
            ux2[i, j] = ux1[i, j] - prm.dt * dpdx
        end
    end

    # add x-direction boundary condition
    boundary_ux(ux2, nx1, nx2, ny1, ny2)

    # the range of y-direction velocity components not including the boundaries
    nx1 = 2
    nx2 = prm.nx + 1
    ny1 = 2
    ny2 = prm.ny + 2

    # y-direction pressureless velocity
    for j in ny1:ny2
        for i in nx1:nx2
            # calculate dpdx
            dpdy = ddx1_1(pp[i, j-1], pp[i, j], dy)
            uy2[i, j] = uy1[i, j] - prm.dt * dpdy
        end
    end

    # add y-direction boundary condition
    boundary_uy(uy2, nx1, nx2, ny1, ny2)

end


end
