using GLMakie, CairoMakie
using DelimitedFiles

filename = "../build/evolution.dqmc.dat"
data = readdlm(filename, comments = true)

file = readlines(filename)
xmin = 0.
xmax = 0.
nbins = 0

for line in file
    if(contains(line, "xmin"))
        xmin = parse(Float64, split(line, ":")[2])
    end
    if(contains(line, "xmax"))
        xmax = parse(Float64, split(line, ":")[2])
    end
    if(contains(line, "nbins"))
        nbins = parse(Int64, split(line, ":")[2])
    end
end

x = y = LinRange(xmin, xmax, nbins)

##

GLMakie.activate!()

with_theme(theme_latexfonts()) do
    fig = Figure();
    ax = Axis3(fig[1,1], 
               xlabel = "x [a.u]",
               ylabel = "y [a.u]")

    hidezdecorations!(ax)

    cm = surface!(ax, x, y, data,
             colormap = :coolwarm,
             label = "Ψ(x,y)")

    ticks_psi = [ round(val, digits = 3) for val in collect(LinRange(findmin(data)[1], findmax(data)[1], 6)) ]

    Colorbar(fig[2,1], cm, label = "Ψ(x,y)", vertical = false, ticks = ticks_psi[2:end-1])

    xlims!(ax, (xmin, xmax))
    ylims!(ax, (xmin, xmax))
    zlims!(ax, (0, nothing))

    display(fig)
end
