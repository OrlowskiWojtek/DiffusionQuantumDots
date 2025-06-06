using GLMakie, CairoMakie
using DelimitedFiles


#filename = "../build/results/2d_harmonic_oscillator/1st_state_2d_anisotropic"
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

#

GLMakie.activate!()

with_theme(theme_latexfonts()) do
    fig = Figure(size = (1024, 768));
    ax = Axis3(fig[1,1], 
               xlabel = "x [nm]",
               ylabel = "y [nm]")
 
    hidezdecorations!(ax)

    cm = surface!(ax, x, y, data,
             colormap = :coolwarm,
             label = "Końcowy rozkład wędrowców")

    Colorbar(fig[2,1], cm, label = "Ψ(x,y)", vertical = false)

    xlims!(ax, (xmin, xmax))
    ylims!(ax, (xmin, xmax))

    display(fig)
    save("plots/1d_2el_excited.png", fig)
end

##

dx = x[2] - x[1]
s = sum(data.^2) * (dx / 0.0529)
