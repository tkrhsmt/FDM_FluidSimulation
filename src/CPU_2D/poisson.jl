# --------------------------------------------------
# MODULE    Poisson
# julia code
#
# DATE : 2024/8/24
# --------------------------------------------------

module Poisson

include("input_param.jl")
using .Param

export poisson!

# 2nd-order central difference
ddx1_1(u1, u2, dx) = (u2 - u1) / (dx)

function poisson!(ux, uy, pp1, pp2, div, prm)

    # the range of pressure components not including the boundaries
    nx1 = 2
    nx2 = prm.nx + 1
    ny1 = 2
    ny2 = prm.ny + 1

    # spatial step size
    dx = prm.lx / prm.nx
    dy = prm.ly / prm.ny

    # calculate velocity divergence
    for j in ny1:ny2
        for i in nx1:nx2
            # du/dx
            tmp1 = ddx1_1(ux[i, j], ux[i+1, j], dx)
            # dv/dy
            tmp2 = ddx1_1(uy[i, j], uy[i, j+1], dy)
            # du/dx + dv/dy
            div[i, j] = tmp1 + tmp2
        end
    end

    # divide by time span
    div /= prm.dt

    # ---------- POISSON METHOD ---------- #

    # convergence conditions for iterative methods
    eps = 1.0e-6
    # max repeat of iterative methods
    maxrep = 1000
    # SOR accelate
    α = 1.4

    for rep in 1:maxrep

        # calculate poisson equation (iterative method) pp1 -> pp2
        SOR_method!(pp1, pp2, div, nx1, nx2, ny1, ny2, dx, dy, α)

        # add boundary
        boundary_pp(pp2, nx1, nx2, ny1, ny2)

        # calculate poisson equation (iterative method) pp2 -> pp1
        SOR_method!(pp2, pp1, div, nx1, nx2, ny1, ny2, dx, dy, α)

        # add boundary
        boundary_pp(pp1, nx1, nx2, ny1, ny2)

        # check residual
        check = check_residual(pp1, div, eps, nx1, nx2, ny1, ny2, dx, dy)

        if check
            break
        end
    end

end

function SOR_method!(pp1, pp2, div, nx1, nx2, ny1, ny2, dx, dy, α)

    # constants used in solving the Poisson equation
    pdxdy = 2.0 * (1 / dx^2 + 1 / dy^2)
    dx2 = 1 / dx^2
    dy2 = 1 / dy^2

    for j in ny1:ny2
        for i in nx1:nx2
            tmp1 = (pp1[i+1, j] + pp1[i-1, j]) * dx2 + (pp1[i, j+1] + pp1[i, j-1]) * dy2
            pp2[i, j] = (tmp1 - div[i, j]) / pdxdy
        end
    end

end

function check_residual(pp, div, eps, nx1, nx2, ny1, ny2, dx, dy)

    # convergence test variables
    ϵ = 0.0

    for j in ny1:ny2
        for i in nx1:nx2
            # calculate residual
            dpdx2 = (pp[i+1, j] - 2 * pp[i, j] + pp[i-1, j]) / (dx^2)
            dpdy2 = (pp[i, j+1] - 2 * pp[i, j] + pp[i, j-1]) / (dy^2)
            res = abs(dpdx2 + dpdy2 - div[i, j])

            # set max residual
            if res > ϵ
                ϵ = res
            end
        end
    end

    # check residual
    if eps > ϵ
        return true
    else
        return false
    end

end

end
