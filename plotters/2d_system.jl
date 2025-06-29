using GLMakie, CairoMakie
using DelimitedFiles

filename = "../data/2el_2d_excited_3mev_5mev_guided/evolution.dqmc.dat"
#filename = "../build/evolution.dqmc.dat"
#filename = "../build/total_distribution.dqmc.dat"

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

CairoMakie.activate!()

with_theme(theme_latexfonts()) do
    fig = Figure(dpi = 300);
    Label(fig[0,1:2], "Końcowy rozkład wędrowców")
    ax = Axis3(fig[1,1], 
               xlabel = L"$x_1$ [nm]",
               ylabel = L"$x_2$ [nm]")

    ax_heatmap = Axis(fig[1,2], aspect = 1,
                      xlabel = L"$x_1$ [nm]",
                      ylabel = L"$x_2$ [nm]")

    hidezdecorations!(ax)

    hm = heatmap!(ax_heatmap, x, y, data,
                  colormap = :coolwarm)

    cm = surface!(ax, x, y, data,
                  colormap = :coolwarm)

    Colorbar(fig[2,1:2], cm, label = L"$f$(\textbf{R})", vertical = false)

    xlims!(ax, (xmin, xmax))
    ylims!(ax, (xmin, xmax))

    xlims!(ax_heatmap, (xmin, xmax))
    ylims!(ax_heatmap, (xmin, xmax))

    display(fig)
    save("plots/2d_2el_excited_3meV_5meV_x1x2.png", fig)
end

##

dx = x[2] - x[1]
s = sum(data.^2) * (dx / 0.0529)^2
