# --------------------------------------------------
# MODULE    les_model
# julia code
#
# DATE : 2024/9/14
# --------------------------------------------------


module Les

export les_model

ddx1_1(u1, u2, dx) = (u2 - u1) / (dx)

function les_model(ux, uy, uz, prm)

    if prm.les_model == 1
        # smagorinsky model
        return smagorinsky_model(ux, uy, uz, prm)
    else prm.les_model == 2
        # coherent structure model
        return coherent_structure_model(ux, uy, uz, prm)
    end

end

function smagorinsky_model(ux, uy, uz, prm)

    dudx = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dudy = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dudz = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dvdx = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dvdy = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dvdz = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dwdx = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dwdy = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dwdz = zeros(prm.nx+2, prm.ny+2, prm.nz+2)

    # spatial step size
    dx = prm.lx / prm.nx
    dy = prm.ly / prm.ny
    dz = prm.lz / prm.nz

    nx1 = 2
    nx2 = prm.nx+1
    ny1 = 2
    ny2 = prm.ny+1
    nz1 = 2
    nz2 = prm.nz+1

    for k in nz1:nz2
        for j in ny1:ny2
            for i in nx1:nx2
                dudx[i, j, k] = ddx1_1(ux[i, j, k], ux[i+1, j, k], dx)
                dudy[i, j, k] = 0.25*(ddx1_1(ux[i, j, k], ux[i, j+1, k], dy) + ddx1_1(ux[i, j-1, k], ux[i, j, k], dy) + ddx1_1(ux[i+1, j, k], ux[i+1, j+1, k], dy) + ddx1_1(ux[i+1, j-1, k], ux[i+1, j, k], dy))
                dudz[i, j, k] = 0.25*(ddx1_1(ux[i, j, k], ux[i, j, k+1], dz) + ddx1_1(ux[i, j, k-1], ux[i, j, k], dz) + ddx1_1(ux[i+1, j, k], ux[i+1, j, k+1], dz) + ddx1_1(ux[i+1, j, k-1], ux[i+1, j, k], dz))
                dvdx[i, j, k] = 0.25*(ddx1_1(uy[i, j, k], uy[i+1, j, k], dx) + ddx1_1(uy[i-1, j, k], uy[i, j, k], dx) + ddx1_1(uy[i, j+1, k], uy[i+1, j+1, k], dx) + ddx1_1(uy[i-1, j+1, k], uy[i, j+1, k], dx))
                dvdy[i, j, k] = ddx1_1(uy[i, j, k], uy[i, j+1, k], dy)
                dvdz[i, j, k] = 0.25*(ddx1_1(uy[i, j, k], uy[i, j, k+1], dz) + ddx1_1(uy[i, j, k-1], uy[i, j, k], dz) + ddx1_1(uy[i, j+1, k], uy[i, j+1, k+1], dz) + ddx1_1(uy[i, j+1, k-1], uy[i, j+1, k], dz))
                dwdx[i, j, k] = 0.25*(ddx1_1(uz[i, j, k], uz[i+1, j, k], dx) + ddx1_1(uz[i-1, j, k], uz[i, j, k], dx) + ddx1_1(uz[i, j, k+1], uz[i+1, j, k+1], dx) + ddx1_1(uz[i-1, j, k+1], uz[i, j, k+1], dx))
                dwdy[i, j, k] = 0.25*(ddx1_1(uz[i, j, k], uz[i, j+1, k], dy) + ddx1_1(uz[i, j-1, k], uz[i, j, k], dy) + ddx1_1(uz[i, j, k+1], uz[i, j+1, k+1], dy) + ddx1_1(uz[i, j-1, k+1], uz[i, j, k+1], dy))
                dwdz[i, j, k] = ddx1_1(uz[i, j, k], uz[i, j, k+1], dz)
            end
        end
    end

    # Deformation Gradient Tensor
    s11 = 0.5 * (dudx + dudx)
    s12 = 0.5 * (dudy + dvdx)
    s13 = 0.5 * (dudz + dwdx)
    s22 = 0.5 * (dvdy + dvdy)
    s23 = 0.5 * (dvdz + dwdy)
    s33 = 0.5 * (dwdz + dwdz)
    ss = s11.*s11 + s22.*s22 + s33.*s33 + 2.0*(s12.*s12 + s13.*s13 + s23.*s23)

    # Filter scale
    Δ = sqrt(dx^2 + dy^2 + dz^2)

    # Smagorinsky constant
    Cs = 0.17

    #　Eddy viscosity coefficient
    ν_t = (Cs*Δ)^2 * sqrt.(2.0*ss)

    # SGS stress divergence
    dτ1 = zeros(prm.nx+3, prm.ny+2, prm.nz+2)
    dτ2 = zeros(prm.nx+2, prm.ny+3, prm.nz+2)
    dτ3 = zeros(prm.nx+2, prm.ny+2, prm.nz+3)

    for k in nz1:nz2
        for j in ny1:ny2
            for i in nx1:nx2+1
                # x-direction SGS stress divergence
                dτ1dx = ddx1_1(-2.0*ν_t[i-1, j, k]*s11[i-1, j, k], -2.0*ν_t[i, j, k]*s11[i, j, k], dx)
                dτ2dy = ddx1_1(-2.0*ν_t[i, j-1, k]*s12[i, j-1, k], -2.0*ν_t[i, j, k]*s12[i, j, k], dy)
                dτ3dz = ddx1_1(-2.0*ν_t[i, j, k-1]*s13[i, j, k-1], -2.0*ν_t[i, j, k]*s13[i, j, k], dz)
                dτ1[i, j, k] = dτ1dx + dτ2dy + dτ3dz
            end
        end
    end
    for k in nz1:nz2
        for j in ny1:ny2+1
            for i in nx1:nx2
                # y-direction SGS stress divergence
                dτ1dx = ddx1_1(-2.0*ν_t[i-1, j, k]*s12[i-1, j, k], -2.0*ν_t[i, j, k]*s12[i, j, k], dx)
                dτ2dy = ddx1_1(-2.0*ν_t[i, j-1, k]*s22[i, j-1, k], -2.0*ν_t[i, j, k]*s22[i, j, k], dy)
                dτ3dz = ddx1_1(-2.0*ν_t[i, j, k-1]*s23[i, j, k-1], -2.0*ν_t[i, j, k]*s23[i, j, k], dz)
                dτ2[i, j, k] = dτ1dx + dτ2dy + dτ3dz
            end
        end
    end
    for k in nz1:nz2+1
        for j in ny1:ny2
            for i in nx1:nx2
                # z-direction SGS stress divergence
                dτ1dx = ddx1_1(-2.0*ν_t[i-1, j, k]*s13[i-1, j, k], -2.0*ν_t[i, j, k]*s13[i, j, k], dx)
                dτ2dy = ddx1_1(-2.0*ν_t[i, j-1, k]*s23[i, j-1, k], -2.0*ν_t[i, j, k]*s23[i, j, k], dy)
                dτ3dz = ddx1_1(-2.0*ν_t[i, j, k-1]*s33[i, j, k-1], -2.0*ν_t[i, j, k]*s33[i, j, k], dz)
                dτ3[i, j, k] = dτ1dx + dτ2dy + dτ3dz
            end
        end
    end

    # x direction boundary
    if prm.bx == 0
        dτ1[1, :, :] = dτ1[end-1, :, :]
        dτ1[end, :, :] = dτ1[2, :, :]
        dτ2[1, :, :] = dτ2[end-1, :, :]
        dτ2[end, :, :] = dτ2[2, :, :]
        dτ3[1, :, :] = dτ3[end-1, :, :]
        dτ3[end, :, :] = dτ3[2, :, :]
    else
        dτ1[1, :, :] = dτ1[2, :, :]
        dτ1[end, :, :] = dτ1[end-1, :, :]
        dτ2[1, :, :] = dτ2[2, :, :]
        dτ2[end, :, :] = dτ2[end-1, :, :]
        dτ3[1, :, :] = dτ3[2, :, :]
        dτ3[end, :, :] = dτ3[end-1, :, :]
    end

    # y direction boundary
    if prm.by == 0
        dτ1[:, 1, :] = dτ1[:, end-1, :]
        dτ1[:, end, :] = dτ1[:, 2, :]
        dτ2[:, 1, :] = dτ2[:, end-1, :]
        dτ2[:, end, :] = dτ2[:, 2, :]
        dτ3[:, 1, :] = dτ3[:, end-1, :]
        dτ3[:, end, :] = dτ3[:, 2, :]
    else
        dτ1[:, 1, :] = dτ1[:, 2, :]
        dτ1[:, end, :] = dτ1[:, end-1, :]
        dτ2[:, 1, :] = dτ2[:, 2, :]
        dτ2[:, end, :] = dτ2[:, end-1, :]
        dτ3[:, 1, :] = dτ3[:, 2, :]
        dτ3[:, end, :] = dτ3[:, end-1, :]
    end

    # z direction boundary
    if prm.bz == 0
        dτ1[:, :, 1] = dτ1[:, :, end-1]
        dτ1[:, :, end] = dτ1[:, :, 2]
        dτ2[:, :, 1] = dτ2[:, :, end-1]
        dτ2[:, :, end] = dτ2[:, :, 2]
        dτ3[:, :, 1] = dτ3[:, :, end-1]
        dτ3[:, :, end] = dτ3[:, :, 2]
    else
        dτ1[:, :, 1] = dτ1[:, :, 2]
        dτ1[:, :, end] = dτ1[:, :, end-1]
        dτ2[:, :, 1] = dτ2[:, :, 2]
        dτ2[:, :, end] = dτ2[:, :, end-1]
        dτ3[:, :, 1] = dτ3[:, :, 2]
        dτ3[:, :, end] = dτ3[:, :, end-1]
    end

    return dτ1, dτ2, dτ3

end

function coherent_structure_model(ux, uy, uz, prm)

    dudx = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dudy = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dudz = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dvdx = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dvdy = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dvdz = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dwdx = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dwdy = zeros(prm.nx+2, prm.ny+2, prm.nz+2)
    dwdz = zeros(prm.nx+2, prm.ny+2, prm.nz+2)

    # spatial step size
    dx = prm.lx / prm.nx
    dy = prm.ly / prm.ny
    dz = prm.lz / prm.nz

    nx1 = 2
    nx2 = prm.nx+1
    ny1 = 2
    ny2 = prm.ny+1
    nz1 = 2
    nz2 = prm.nz+1

    for k in nz1:nz2
        for j in ny1:ny2
            for i in nx1:nx2
                dudx[i, j, k] = ddx1_1(ux[i, j, k], ux[i+1, j, k], dx)
                dudy[i, j, k] = 0.25*(ddx1_1(ux[i, j, k], ux[i, j+1, k], dy) + ddx1_1(ux[i, j-1, k], ux[i, j, k], dy) + ddx1_1(ux[i+1, j, k], ux[i+1, j+1, k], dy) + ddx1_1(ux[i+1, j-1, k], ux[i+1, j, k], dy))
                dudz[i, j, k] = 0.25*(ddx1_1(ux[i, j, k], ux[i, j, k+1], dz) + ddx1_1(ux[i, j, k-1], ux[i, j, k], dz) + ddx1_1(ux[i+1, j, k], ux[i+1, j, k+1], dz) + ddx1_1(ux[i+1, j, k-1], ux[i+1, j, k], dz))
                dvdx[i, j, k] = 0.25*(ddx1_1(uy[i, j, k], uy[i+1, j, k], dx) + ddx1_1(uy[i-1, j, k], uy[i, j, k], dx) + ddx1_1(uy[i, j+1, k], uy[i+1, j+1, k], dx) + ddx1_1(uy[i-1, j+1, k], uy[i, j+1, k], dx))
                dvdy[i, j, k] = ddx1_1(uy[i, j, k], uy[i, j+1, k], dy)
                dvdz[i, j, k] = 0.25*(ddx1_1(uy[i, j, k], uy[i, j, k+1], dz) + ddx1_1(uy[i, j, k-1], uy[i, j, k], dz) + ddx1_1(uy[i, j+1, k], uy[i, j+1, k+1], dz) + ddx1_1(uy[i, j+1, k-1], uy[i, j+1, k], dz))
                dwdx[i, j, k] = 0.25*(ddx1_1(uz[i, j, k], uz[i+1, j, k], dx) + ddx1_1(uz[i-1, j, k], uz[i, j, k], dx) + ddx1_1(uz[i, j, k+1], uz[i+1, j, k+1], dx) + ddx1_1(uz[i-1, j, k+1], uz[i, j, k+1], dx))
                dwdy[i, j, k] = 0.25*(ddx1_1(uz[i, j, k], uz[i, j+1, k], dy) + ddx1_1(uz[i, j-1, k], uz[i, j, k], dy) + ddx1_1(uz[i, j, k+1], uz[i, j+1, k+1], dy) + ddx1_1(uz[i, j-1, k+1], uz[i, j, k+1], dy))
                dwdz[i, j, k] = ddx1_1(uz[i, j, k], uz[i, j, k+1], dz)
            end
        end
    end

    # Deformation Gradient Tensor
    s11 = 0.5 * (dudx + dudx)
    s12 = 0.5 * (dudy + dvdx)
    s13 = 0.5 * (dudz + dwdx)
    s22 = 0.5 * (dvdy + dvdy)
    s23 = 0.5 * (dvdz + dwdy)
    s33 = 0.5 * (dwdz + dwdz)
    ss = s11.*s11 + s22.*s22 + s33.*s33 + 2.0*(s12.*s12 + s13.*s13 + s23.*s23)

    # Rotation Tensor
    w12 = 0.5 * (dudy - dvdx)
    w13 = 0.5 * (dudz - dwdx)
    w23 = 0.5 * (dvdz - dwdy)
    ww = 2.0*(w12.*w12 + w13.*w13 + w23.*w23)

    # Filter scale
    Δ = sqrt(dx^2 + dy^2 + dz^2)

    # coherent structure function
    Q = 0.5 * (ww - ss)
    E = 0.5 * (ww + ss)
    Fcs = Q ./ (E .+ 1.0e-6)

    # Smagorinsky constant
    C = 1/20 * abs.(Fcs).^(3/2)

    #　Eddy viscosity coefficient
    ν_t = Δ^2 * C .* sqrt.(2.0*ss)

    # SGS stress divergence
    dτ1 = zeros(prm.nx+3, prm.ny+2, prm.nz+2)
    dτ2 = zeros(prm.nx+2, prm.ny+3, prm.nz+2)
    dτ3 = zeros(prm.nx+2, prm.ny+2, prm.nz+3)

    for k in nz1:nz2
        for j in ny1:ny2
            for i in nx1:nx2+1
                # x-direction SGS stress divergence
                dτ1dx = ddx1_1(-2.0*ν_t[i-1, j, k]*s11[i-1, j, k], -2.0*ν_t[i, j, k]*s11[i, j, k], dx)
                dτ2dy = ddx1_1(-2.0*ν_t[i, j-1, k]*s12[i, j-1, k], -2.0*ν_t[i, j, k]*s12[i, j, k], dy)
                dτ3dz = ddx1_1(-2.0*ν_t[i, j, k-1]*s13[i, j, k-1], -2.0*ν_t[i, j, k]*s13[i, j, k], dz)
                dτ1[i, j, k] = dτ1dx + dτ2dy + dτ3dz
            end
        end
    end
    for k in nz1:nz2
        for j in ny1:ny2+1
            for i in nx1:nx2
                # y-direction SGS stress divergence
                dτ1dx = ddx1_1(-2.0*ν_t[i-1, j, k]*s12[i-1, j, k], -2.0*ν_t[i, j, k]*s12[i, j, k], dx)
                dτ2dy = ddx1_1(-2.0*ν_t[i, j-1, k]*s22[i, j-1, k], -2.0*ν_t[i, j, k]*s22[i, j, k], dy)
                dτ3dz = ddx1_1(-2.0*ν_t[i, j, k-1]*s23[i, j, k-1], -2.0*ν_t[i, j, k]*s23[i, j, k], dz)
                dτ2[i, j, k] = dτ1dx + dτ2dy + dτ3dz
            end
        end
    end
    for k in nz1:nz2+1
        for j in ny1:ny2
            for i in nx1:nx2
                # z-direction SGS stress divergence
                dτ1dx = ddx1_1(-2.0*ν_t[i-1, j, k]*s13[i-1, j, k], -2.0*ν_t[i, j, k]*s13[i, j, k], dx)
                dτ2dy = ddx1_1(-2.0*ν_t[i, j-1, k]*s23[i, j-1, k], -2.0*ν_t[i, j, k]*s23[i, j, k], dy)
                dτ3dz = ddx1_1(-2.0*ν_t[i, j, k-1]*s33[i, j, k-1], -2.0*ν_t[i, j, k]*s33[i, j, k], dz)
                dτ3[i, j, k] = dτ1dx + dτ2dy + dτ3dz
            end
        end
    end

    # x direction boundary
    if prm.bx == 0
        dτ1[1, :, :] = dτ1[end-1, :, :]
        dτ1[end, :, :] = dτ1[2, :, :]
        dτ2[1, :, :] = dτ2[end-1, :, :]
        dτ2[end, :, :] = dτ2[2, :, :]
        dτ3[1, :, :] = dτ3[end-1, :, :]
        dτ3[end, :, :] = dτ3[2, :, :]
    else
        dτ1[1, :, :] = dτ1[2, :, :]
        dτ1[end, :, :] = dτ1[end-1, :, :]
        dτ2[1, :, :] = dτ2[2, :, :]
        dτ2[end, :, :] = dτ2[end-1, :, :]
        dτ3[1, :, :] = dτ3[2, :, :]
        dτ3[end, :, :] = dτ3[end-1, :, :]
    end

    # y direction boundary
    if prm.by == 0
        dτ1[:, 1, :] = dτ1[:, end-1, :]
        dτ1[:, end, :] = dτ1[:, 2, :]
        dτ2[:, 1, :] = dτ2[:, end-1, :]
        dτ2[:, end, :] = dτ2[:, 2, :]
        dτ3[:, 1, :] = dτ3[:, end-1, :]
        dτ3[:, end, :] = dτ3[:, 2, :]
    else
        dτ1[:, 1, :] = dτ1[:, 2, :]
        dτ1[:, end, :] = dτ1[:, end-1, :]
        dτ2[:, 1, :] = dτ2[:, 2, :]
        dτ2[:, end, :] = dτ2[:, end-1, :]
        dτ3[:, 1, :] = dτ3[:, 2, :]
        dτ3[:, end, :] = dτ3[:, end-1, :]
    end

    # z direction boundary
    if prm.bz == 0
        dτ1[:, :, 1] = dτ1[:, :, end-1]
        dτ1[:, :, end] = dτ1[:, :, 2]
        dτ2[:, :, 1] = dτ2[:, :, end-1]
        dτ2[:, :, end] = dτ2[:, :, 2]
        dτ3[:, :, 1] = dτ3[:, :, end-1]
        dτ3[:, :, end] = dτ3[:, :, 2]
    else
        dτ1[:, :, 1] = dτ1[:, :, 2]
        dτ1[:, :, end] = dτ1[:, :, end-1]
        dτ2[:, :, 1] = dτ2[:, :, 2]
        dτ2[:, :, end] = dτ2[:, :, end-1]
        dτ3[:, :, 1] = dτ3[:, :, 2]
        dτ3[:, :, end] = dτ3[:, :, end-1]
    end

    return dτ1, dτ2, dτ3

end

end
