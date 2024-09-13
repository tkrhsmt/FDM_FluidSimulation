# --------------------------------------------------
# MODULE    Poisson
# julia code
#
# DATE : 2024/8/24
# --------------------------------------------------

module Poisson

using FFTW

export poisson!, poisson_fft

# 2nd-order central difference
ddx1_1(u1, u2, dx) = (u2 - u1) / (dx)

function poisson_fft(ux, uy, uz, pp1, pp2, divu, prm)

    # the range of pressure components not including the boundaries
    nx1 = 2
    nx2 = prm.nx + 1
    ny1 = 2
    ny2 = prm.ny + 1
    nz1 = 2
    nz2 = prm.nz + 1

    # spatial step size
    dx = prm.lx / prm.nx
    dy = prm.ly / prm.ny
    dz = prm.lz / prm.nz

    # calculate velocity divergence
    for k in nz1:nz2
        for j in ny1:ny2
            for i in nx1:nx2
                # du/dx
                tmp1 = ddx1_1(ux[i, j, k], ux[i+1, j, k], dx)
                # dv/dy
                tmp2 = ddx1_1(uy[i, j, k], uy[i, j+1, k], dy)
                # dw/dz
                tmp3 = ddx1_1(uz[i, j, k], uz[i, j, k+1], dz)
                # du/dx + dv/dy + dw/dz
                divu[i, j, k] = tmp1 + tmp2 + tmp3
            end
        end
    end

    # divide by time span
    divu /= prm.dt

    # wave number (periodic boundary)
    kx = fftfreq(prm.nx)
    ky = fftfreq(prm.ny)
    kz = fftfreq(prm.nz)


    # calculate pp_f
    if prm.bx == 0 && prm.by == 0 && prm.bz == 0

        pp_f = zeros(Complex{Float64}, prm.nx, prm.ny, prm.nz)

        divu_f = FFTW.fft(divu[2:end-1, 2:end-1, 2:end-1])
        pp_f = fft_calculate_000!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)
        pp1[2:end-1, 2:end-1, 2:end-1] = real(FFTW.ifft(pp_f))
        boundary_pp_000!(pp1)

    elseif prm.bx == 0 && prm.by == 0 && prm.bz == 1

        pp_f = zeros(Complex{Float64}, prm.nx, prm.ny, prm.nz)

        divu_f1 = FFTW.r2r(divu[2:end-1, 2:end-1, 2:end-1],FFTW.REDFT00, 3)
        divu_f1 = FFTW.fft(divu_f1, 2)
        divu_f = FFTW.fft(divu_f1, 1)

        pp_f = fft_calculate_001!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)

        pp_f1 = FFTW.ifft(pp_f, 1)
        pp_f1 = real(FFTW.ifft(pp_f1, 2))
        pp1[2:end-1, 2:end-1, 2:end-1] = FFTW.r2r(pp_f1,FFTW.REDFT00, 3)/(2*prm.nz)
        boundary_pp_001!(pp1)

    elseif prm.bx == 0 && prm.by == 1 && prm.bz == 0

        pp_f = zeros(Complex{Float64}, prm.nx, prm.ny, prm.nz)

        divu_f1 = FFTW.r2r(divu[2:end-1, 2:end-1, 2:end-1],FFTW.REDFT00, 2)
        divu_f1 = FFTW.fft(divu_f1, 3)
        divu_f = FFTW.fft(divu_f1, 1)

        pp_f = fft_calculate_010!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)

        pp_f1 = FFTW.ifft(pp_f, 1)
        pp_f1 = real(FFTW.ifft(pp_f1, 3))
        pp1[2:end-1, 2:end-1, 2:end-1] = FFTW.r2r(pp_f1,FFTW.REDFT00, 2)/(2*prm.ny)
        boundary_pp_010!(pp1)

    elseif prm.bx == 1 && prm.by == 0 && prm.bz == 0

        pp_f = zeros(Complex{Float64}, prm.nx, prm.ny, prm.nz)

        divu_f1 = FFTW.r2r(divu[2:end-1, 2:end-1, 2:end-1],FFTW.REDFT00, 1)
        divu_f1 = FFTW.fft(divu_f1, 2)
        divu_f = FFTW.fft(divu_f1, 3)

        pp_f = fft_calculate_100!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)

        pp_f1 = FFTW.ifft(pp_f, 3)
        pp_f1 = real(FFTW.ifft(pp_f1, 2))
        pp1[2:end-1, 2:end-1, 2:end-1] = FFTW.r2r(pp_f1,FFTW.REDFT00, 1)/(2*prm.nx)
        boundary_pp_100!(pp1)

    elseif prm.bx == 1 && prm.by == 1 && prm.bz == 0

        pp_f = zeros(Complex{Float64}, prm.nx, prm.ny, prm.nz)

        divu_f1 = FFTW.r2r(divu[2:end-1, 2:end-1, 2:end-1],FFTW.REDFT00, 1)
        divu_f1 = FFTW.r2r(divu_f1,FFTW.REDFT00, 2)
        divu_f = FFTW.fft(divu_f1, 3)

        pp_f = fft_calculate_110!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)

        pp_f1 = real(FFTW.ifft(pp_f, 3))
        pp_f1 = FFTW.r2r(pp_f1,FFTW.REDFT00, 2)/(2*prm.ny)
        pp1[2:end-1, 2:end-1, 2:end-1] = FFTW.r2r(pp_f1,FFTW.REDFT00, 1)/(2*prm.nx)
        boundary_pp_110!(pp1)

    elseif prm.bx == 1 && prm.by == 0 && prm.bz == 1

        pp_f = zeros(Complex{Float64}, prm.nx, prm.ny, prm.nz)

        divu_f1 = FFTW.r2r(divu[2:end-1, 2:end-1, 2:end-1],FFTW.REDFT00, 1)
        divu_f1 = FFTW.r2r(divu_f1,FFTW.REDFT00, 3)
        divu_f = FFTW.fft(divu_f1, 2)

        pp_f = fft_calculate_101!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)

        pp_f1 = real(FFTW.ifft(pp_f, 2))
        pp_f1 = FFTW.r2r(pp_f1,FFTW.REDFT00, 3)/(2*prm.nz)
        pp1[2:end-1, 2:end-1, 2:end-1] = FFTW.r2r(pp_f1,FFTW.REDFT00, 1)/(2*prm.nx)
        boundary_pp_101!(pp1)

    elseif prm.bx == 0 && prm.by == 1 && prm.bz == 1

        pp_f = zeros(Complex{Float64}, prm.nx, prm.ny, prm.nz)

        divu_f1 = FFTW.r2r(divu[2:end-1, 2:end-1, 2:end-1],FFTW.REDFT00, 3)
        divu_f1 = FFTW.r2r(divu_f1,FFTW.REDFT00, 2)
        divu_f = FFTW.fft(divu_f1, 1)

        pp_f = fft_calculate_011!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)

        pp_f1 = real(FFTW.ifft(pp_f, 1))
        pp_f1 = FFTW.r2r(pp_f1,FFTW.REDFT00, 2)/(2*prm.ny)
        pp1[2:end-1, 2:end-1, 2:end-1] = FFTW.r2r(pp_f1,FFTW.REDFT00, 3)/(2*prm.nz)
        boundary_pp_011!(pp1)

    elseif prm.bx == 1 && prm.by == 1 && prm.bz == 1

        pp_f = zeros(Complex{Float64}, prm.nx, prm.ny, prm.nz)

        divu_f = FFTW.r2r(divu[2:end-1, 2:end-1, 2:end-1],FFTW.REDFT00)

        pp_f = fft_calculate_111!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)

        pp1[2:end-1, 2:end-1, 2:end-1] = FFTW.r2r(pp_f,FFTW.REDFT00)/((2*prm.nx)*(2*prm.ny)*(2*prm.nz))
        boundary_pp_111!(pp1)
    end

    return 1

end

function boundary_pp_000!(pp)

    pp[1, :, :] = pp[end-1, :, :]
    pp[end, :, :] = pp[2, :, :]
    pp[:, 1, :] = pp[:, end-1, :]
    pp[:, end, :] = pp[:, 2, :]
    pp[:, :, 1] = pp[:, :, end-1]
    pp[:, :, end] = pp[:, :, 2]

end

function boundary_pp_100!(pp)

    pp[1, :, :] = pp[2, :, :]
    pp[end, :, :] = pp[end-1, :, :]
    pp[:, 1, :] = pp[:, end-1, :]
    pp[:, end, :] = pp[:, 2, :]
    pp[:, :, 1] = pp[:, :, end-1]
    pp[:, :, end] = pp[:, :, 2]

end

function boundary_pp_010!(pp)

    pp[1, :, :] = pp[end-1, :, :]
    pp[end, :, :] = pp[2, :, :]
    pp[:, 1, :] = pp[:, 2, :]
    pp[:, end, :] = pp[:, end-1, :]
    pp[:, :, 1] = pp[:, :, end-1]
    pp[:, :, end] = pp[:, :, 2]

end

function boundary_pp_001!(pp)

    pp[1, :, :] = pp[end-1, :, :]
    pp[end, :, :] = pp[2, :, :]
    pp[:, 1, :] = pp[:, end-1, :]
    pp[:, end, :] = pp[:, 2, :]
    pp[:, :, 1] = pp[:, :, 2]
    pp[:, :, end] = pp[:, :, end-1]

end

function boundary_pp_110!(pp)

    pp[1, :, :] = pp[2, :, :]
    pp[end, :, :] = pp[end-1, :, :]
    pp[:, 1, :] = pp[:, 2, :]
    pp[:, end, :] = pp[:, end-1, :]
    pp[:, :, 1] = pp[:, :, end-1]
    pp[:, :, end] = pp[:, :, 2]

end

function boundary_pp_101!(pp)

    pp[1, :, :] = pp[2, :, :]
    pp[end, :, :] = pp[end-1, :, :]
    pp[:, 1, :] = pp[:, end-1, :]
    pp[:, end, :] = pp[:, 2, :]
    pp[:, :, 1] = pp[:, :, 2]
    pp[:, :, end] = pp[:, :, end-1]

end

function boundary_pp_011!(pp)

    pp[1, :, :] = pp[end-1, :, :]
    pp[end, :, :] = pp[2, :, :]
    pp[:, 1, :] = pp[:, 2, :]
    pp[:, end, :] = pp[:, end-1, :]
    pp[:, :, 1] = pp[:, :, 2]
    pp[:, :, end] = pp[:, :, end-1]

end

function boundary_pp_111!(pp)

    pp[1, :, :] = pp[2, :, :]
    pp[end, :, :] = pp[end-1, :, :]
    pp[:, 1, :] = pp[:, 2, :]
    pp[:, end, :] = pp[:, end-1, :]
    pp[:, :, 1] = pp[:, :, 2]
    pp[:, :, end] = pp[:, :, end-1]

end

function fft_calculate_000!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)
    for k in 1:prm.nz
        for j in 1:prm.ny
            for i in 1:prm.nx
                tmp = (2.0/(dx*dx))*(cos(2*pi*kx[i]) - 1.0) + (2.0/(dy*dy))*(cos(2*pi*ky[j]) - 1.0) + (2.0/(dz*dz))*(cos(2*pi*kz[k]) - 1.0)
                pp_f[i, j, k] = divu_f[i, j, k]/(tmp+1.0e-8)
            end
        end
    end
    pp_f[1, 1, 1] = 0.0
    return pp_f
end

function fft_calculate_100!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)
    for k in 1:prm.nz
        for j in 1:prm.ny
            for i in 1:prm.nx
                tmp = (2.0/(dx*dx))*(cos(pi*i/prm.nx) - 1.0) + (2.0/(dy*dy))*(cos(2*pi*ky[j]) - 1.0) + (2.0/(dz*dz))*(cos(2*pi*kz[k]) - 1.0)
                pp_f[i, j, k] = divu_f[i, j, k]/tmp
            end
        end
    end
    return pp_f
end

function fft_calculate_010!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)
    for k in 1:prm.nz
        for j in 1:prm.ny
            for i in 1:prm.nx
                tmp = (2.0/(dx*dx))*(cos(2*pi*kx[i]) - 1.0) + (2.0/(dy*dy))*(cos(pi*j/prm.ny) - 1.0) + (2.0/(dz*dz))*(cos(2*pi*kz[k]) - 1.0)
                pp_f[i, j, k] = divu_f[i, j, k]/tmp
            end
        end
    end
    return pp_f
end

function fft_calculate_001!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)
    for k in 1:prm.nz
        for j in 1:prm.ny
            for i in 1:prm.nx
                tmp = (2.0/(dx*dx))*(cos(2*pi*kx[i]) - 1.0) + (2.0/(dy*dy))*(cos(2*pi*ky[j]) - 1.0) + (2.0/(dz*dz))*(cos(pi*k/prm.nz) - 1.0)
                pp_f[i, j, k] = divu_f[i, j, k]/tmp
            end
        end
    end
    return pp_f
end

function fft_calculate_110!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)
    for k in 1:prm.nz
        for j in 1:prm.ny
            for i in 1:prm.nx
                tmp = (2.0/(dx*dx))*(cos(pi*i/prm.nx) - 1.0) + (2.0/(dy*dy))*(cos(pi*j/prm.ny) - 1.0) + (2.0/(dz*dz))*(cos(2*pi*kz[k]) - 1.0)
                pp_f[i, j, k] = divu_f[i, j, k]/tmp
            end
        end
    end
    return pp_f
end

function fft_calculate_101!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)
    for k in 1:prm.nz
        for j in 1:prm.ny
            for i in 1:prm.nx
                tmp = (2.0/(dx*dx))*(cos(pi*i/prm.nx) - 1.0) + (2.0/(dy*dy))*(cos(2*pi*ky[j]) - 1.0) + (2.0/(dz*dz))*(cos(pi*k/prm.nz) - 1.0)
                pp_f[i, j, k] = divu_f[i, j, k]/tmp
            end
        end
    end
    return pp_f
end

function fft_calculate_011!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)
    for k in 1:prm.nz
        for j in 1:prm.ny
            for i in 1:prm.nx
                tmp = (2.0/(dx*dx))*(cos(2*pi*kx[i]) - 1.0) + (2.0/(dy*dy))*(cos(pi*j/prm.ny) - 1.0) + (2.0/(dz*dz))*(cos(pi*k/prm.nz) - 1.0)
                pp_f[i, j, k] = divu_f[i, j, k]/tmp
            end
        end
    end
    return pp_f
end

function fft_calculate_111!(pp_f, divu_f, prm, dx, dy, dz, kx, ky, kz)
    for k in 1:prm.nz
        for j in 1:prm.ny
            for i in 1:prm.nx
                tmp = (2.0/(dx*dx))*(cos(pi*i/prm.nx) - 1.0) + (2.0/(dy*dy))*(cos(pi*j/prm.ny) - 1.0) + (2.0/(dz*dz))*(cos(pi*k/prm.nz) - 1.0)
                pp_f[i, j, k] = divu_f[i, j, k]/tmp
            end
        end
    end
    return pp_f
end

end
