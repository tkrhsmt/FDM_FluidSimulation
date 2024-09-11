# --------------------------------------------------
# FluidSimulation_CPU constant file
# julia code
# --------------------------------------------------

# number of lattices
NX = 128
NY = 64
NZ = 64

# length of flow field
LX = 2.0
LY = 1.0
LZ = 1.0

# start and end number
ISTART = 1
IEND = 30000

# time span
DT = 0.001

# output log file
ILOG = 10

# output vtk file
IOUTPUT = 100

# reynolds number
RE = 2000

# poisson solver (0 : FFT spectral method, 1 : multigrid method)
POISSON_SCHEME = 0

# pressure boundary (0 : periodic, 1 : neumann)
BX = 1
BY = 0
BZ = 0

# boundary condition ------------------------------

function boundary_ux(ux, nx1, nx2, ny1, ny2, nz1, nz2)

    ux[1:2, :, :] .= 1.0
    ux[end, :, :] = ux[end-1, :, :]
    ux[:, 1, :] = ux[:, end-1, :]
    ux[:, end, :] = ux[:, 2, :]
    ux[:, :, 1] = ux[:, :, end-1]
    ux[:, :, end] = ux[:, :, 2]

    nx = nx2 - nx1 + 1
    ny = ny2 - ny1 + 1
    dx = 2.0 / nx
    dy = 1.0 / ny

    for k in nz1:nz2
        for j in ny1:ny2
            y = dy * (j - ny1) - 0.5
            for i in nx1:nx2
                x = dx * (i - nx1) - 0.5
                l = sqrt(x^2 + y^2)
                if l < 0.1
                    ux[i, j, k] = 0.0
                end
            end
        end
    end

end

function boundary_uy(uy, nx1, nx2, ny1, ny2, nz1, nz2)

    uy[:, 1, :] = uy[:, end-2, :]
    uy[:, end, :] = uy[:, 3, :]
    uy[1, :, :] =  - uy[2, :, :]
    uy[end, :, :] = uy[end-1, :, :]
    uy[:, :, 1] = uy[:, :, end-1]
    uy[:, :, end] = uy[:, :, 2]

    nx = nx2 - nx1 + 1
    ny = ny2 - ny1 + 1
    dx = 2.0 / nx
    dy = 1.0 / ny

    for k in nz1:nz2
        for j in ny1:ny2
            y = dy * (j - ny1) - 0.5
            for i in nx1:nx2
                x = dx * (i - nx1) - 0.5
                l = sqrt(x^2 + y^2)
                if l < 0.1
                    uy[i, j, k] = 0.0
                end
            end
        end
    end

end

function boundary_uz(uz, nx1, nx2, ny1, ny2, nz1, nz2)

    uz[1, :, :] = - uz[2, :, :]
    uz[end, :, :] = uz[end-1, :, :]
    uz[:, 1, :] = uz[:, end-1, :]
    uz[:, end, :] = uz[:, 2, :]
    uz[:, :, 1] = uz[:, :, end-1]
    uz[:, :, end] = uz[:, :, 2]

    nx = nx2 - nx1 + 1
    ny = ny2 - ny1 + 1
    dx = 2.0 / nx
    dy = 1.0 / ny

    for k in nz1:nz2
        for j in ny1:ny2
            y = dy * (j - ny1) - 0.5
            for i in nx1:nx2
                x = dx * (i - nx1) - 0.5
                l = sqrt(x^2 + y^2)
                if l < 0.1
                    uz[i, j, k] = 0.0
                end
            end
        end
    end

end

# init condition ------------------------------

function init!(ux1, uy1, uz1, prm)

    ux1[:,:,:] = 0.001 * rand(prm.nx+3, prm.ny+2, prm.nz+2)
    uy1[:,:,:] = 0.001 * rand(prm.nx+2, prm.ny+3, prm.nz+2)
    uz1[:,:,:] = 0.001 * rand(prm.nx+2, prm.ny+2, prm.nz+3)
end

# velocity force ------------------------------

function force_ux(ux, nx1, nx2, ny1, ny2, nz1, nz2, dt)

end

function force_uy(uy, nx1, nx2, ny1, ny2, nz1, nz2, dt)

end

function force_uz(uz, nx1, nx2, ny1, ny2, nz1, nz2, dt)

end
