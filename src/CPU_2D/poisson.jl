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
    eps = 1.0e-4
    # max repeat of iterative methods
    maxrep = 1000
    # SOR accelate
    α = 1.0
    # repeat counter
    rep_counter = 0

    for rep in 1:maxrep
        rep_counter += 1

        #iterate_method!(pp1, pp2, div, nx1, nx2, ny1, ny2, dx, dy, α)
        multigrid_v_cycle!(pp1, pp2, div, nx1, nx2, ny1, ny2, dx, dy, α; num_cycles=1)

        # check residual
        check = check_residual(pp1, div, eps, nx1, nx2, ny1, ny2, dx, dy)

        if check
            break
        end
    end

    return rep_counter

end

function multigrid_v_cycle!(pp1, pp2, f, nx1, nx2, ny1, ny2, dx, dy, α; num_cycles=1)

    # x-direction spatial size
    nx = nx2 - nx1 + 1
    # y-direction spatial size
    ny = ny2 - ny1 + 1

    for time in 1:num_cycles

        # Pre-smoothing
        iterate_method!(pp1, pp2, f, nx1, nx2, ny1, ny2, dx, dy, α)

        # Compute residual
        res = copy(f)
        compute_grid_residual!(res, pp1, f, nx1, nx2, ny1, ny2, dx, dy)

        # Restrict residual
        res_c, nxc1, nxc2, nyc1, nyc2 = restrict(res, nx, ny)
        u_c1 = zeros(size(res_c))
        u_c2 = zeros(size(res_c))

        # Solve on coarse grid (recursively or directly)
        if nx > 2 && ny > 2
            multigrid_v_cycle!(u_c1, u_c2, res_c, nxc1, nxc2, nyc1, nyc2, 2 * dx, 2 * dy, α, num_cycles=1)
        end

        # Prolongate and correct
        pp1 .+= prolongate(u_c1, nx, ny, nxc1, nxc2, nyc1, nyc2)

        iterate_method!(pp1, pp2, f, nx1, nx2, ny1, ny2, dx, dy, α)

    end
    return pp1
end

function prolongate(uc, nx, ny, nxc1, nxc2, nyc1, nyc2)

    # fine grid error
    u = zeros(nx + 2, ny + 2)

    for i in nxc1:nxc2
        for j in nyc1:nyc2
            # center of grid
            u[2*i-2, 2*j-2] = uc[i, j]
            # edge of grid
            u[2*i+1-2, 2*j-2] = 0.5 * (uc[i, j] + uc[i+1, j])
            u[2*i-2, 2*j+1-2] = 0.5 * (uc[i, j] + uc[i, j+1])
            # corner of grid
            u[2*i+1-2, 2*j+1-2] = 0.25 * (uc[i, j] + uc[i+1, j] + uc[i, j+1] + uc[i+1, j+1])
        end
    end
    return u
end

function restrict(u, nx, ny)

    # rough grid size
    nxc, nyc = div(nx + 1, 2) - 1, div(ny + 1, 2) - 1
    # rough grid residual
    uc = zeros(nxc + 2, nyc + 2)

    # the range of rough grid residual components not including the boundaries
    nxc1 = 2
    nxc2 = nxc + 1
    nyc1 = 2
    nyc2 = nyc + 1

    for i in nxc1:nxc2
        for j in nyc1:nyc2

            # o - * - o
            # |   |   |
            # * - + - *
            # |   |   |
            # o - * - o
            #
            # + : center of grid -> tmp1
            # * : edge of grid   -> tmp2
            # o : corner of grid -> tmp3

            tmp1 = u[2*i-2, 2*j-2]
            tmp2 = u[2*i+1-2, 2*j-2] + u[2*i-1-2, 2*j-2] + u[2*i-2, 2*j+1-2] + u[2*i-2, 2*j-1-2]
            tmp3 = u[2*i+1-2, 2*j+1-2] + u[2*i-1-2, 2*j+1-2] + u[2*i+1-2, 2*j-1-2] + u[2*i-1-2, 2*j-1-2]
            uc[i, j] = (4 * tmp1 + 2 * tmp2 + tmp3) / 16.0
        end
    end
    return uc, nxc1, nxc2, nyc1, nyc2
end

function compute_grid_residual!(res, pp, divx, nx1, nx2, ny1, ny2, dx, dy)

    for j in ny1:ny2
        for i in nx1:nx2
            # calculate residual
            # d/dx(dp/dx)
            dpdx2 = (pp[i+1, j] - 2 * pp[i, j] + pp[i-1, j]) / (dx^2)
            # d/dy(dp/dy)
            dpdy2 = (pp[i, j+1] - 2 * pp[i, j] + pp[i, j-1]) / (dy^2)
            # div(u)/dt - d/dx(dp/dx) - d/dy(dp/dy)
            res[i, j] = divx[i, j] - dpdx2 - dpdy2
        end
    end

    #return maximum(abs.(res))
end

function iterate_method!(pp1, pp2, div, nx1, nx2, ny1, ny2, dx, dy, α)

    for _ in 1:2
        # calculate poisson equation (iterative method) pp1 -> pp2
        SOR_method!(pp1, pp2, div, nx1, nx2, ny1, ny2, dx, dy, α)

        # add boundary
        boundary_pp(pp2, nx1, nx2, ny1, ny2)

        # calculate poisson equation (iterative method) pp2 -> pp1
        SOR_method!(pp2, pp1, div, nx1, nx2, ny1, ny2, dx, dy, α)

        # add boundary
        boundary_pp(pp1, nx1, nx2, ny1, ny2)
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
            pp2[i, j] = α * (tmp1 - div[i, j]) / pdxdy + (1 - α) * pp1[i, j]
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
