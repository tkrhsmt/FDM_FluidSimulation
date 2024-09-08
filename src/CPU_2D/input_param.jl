# --------------------------------------------------
# MODULE    Param
# julia code
#
# DATE : 2024/8/23
# --------------------------------------------------

module Param

using WriteVTK
using HDF5

export inputfile_param, boundary_ux, boundary_uy, boundary_uz, boundary_pp, init!, output_vtkfile, inputdata

struct parameter

    nx :: Int64
    ny :: Int64

    lx :: Float64
    ly :: Float64

    istart :: Int64
    iend :: Int64

    dt :: Float64

    ilog :: Int64
    ioutput :: Int64

    re :: Float64

    bx :: Int64
    by :: Int64

end

function boundary_ux(ux, nx1, nx2, ny1, ny2)

    ux[1:2, :] .= 0.0
    ux[end-1:end, :] .= 0.0
    ux[:, 1] = - ux[:, 2] .+ 2.0
    ux[:, end] = - ux[:, end-1]


end

function boundary_uy(uy, nx1, nx2, ny1, ny2)

    uy[:, 1:2] .= 0.0
    uy[:, end-1:end] .= 0.0
    uy[1, :] = - uy[2, :]
    uy[end, :] = - uy[end-1, :]

end

function boundary_pp(pp, nx1, nx2, ny1, ny2)

    pp[1, :] = pp[2, :]
    pp[end, :] = pp[end-1, :]
    pp[:, 1] = pp[:, 2]
    #pp[:, 1] .= 0.0
    pp[:, end] = pp[:, end-1]

    pav = sum(pp)
    pp = pp .- pav

end

function init!(ux1, uy1, prm)

    ux1 .= 1.0
end

function output_vtkfile(ux, uy, pp, fx, fy, prm, filename)

    dx = prm.lx / prm.nx
    dy = prm.ly / prm.ny

    x = dx:dx:dx*prm.nx
    y = dy:dy:dy*prm.ny

    # output vtk file
    vtk_grid(filename, x, y) do vtk
        vtk["ux"] = 0.5*(ux[2:prm.nx+1, 2:prm.ny+1] + ux[3:prm.nx+2, 2:prm.ny+1])
        vtk["uy"] = 0.5*(uy[2:prm.nx+1, 2:prm.ny+1] + uy[2:prm.nx+1, 3:prm.ny+2])
        vtk["pp"] = pp[2:prm.nx+1, 2:prm.ny+1]
    end

    # output checkout file
    @time h5open("checkpoint", "w") do file
        write(file,"ux",ux)
        write(file,"uy",uy)
        write(file,"fx",fx)
        write(file,"fy",fy)
    end

end

function inputfile_param(INPUT_FILE)

    # include parameter file
    include(INPUT_FILE)
    # setting parameter struct
    prm = parameter(NX, NY, LX, LY, ISTART, IEND, DT, ILOG, IOUTPUT, RE, BX, BY)

    return prm
end

function inputdata()

    file = h5open("checkpoint", "r")
    ux = read(file,"ux")
    uy = read(file,"uy")
    fx = read(file,"fx")
    fy = read(file,"fy")

    close(file)
    return ux, uy, fx, fy
end

end
