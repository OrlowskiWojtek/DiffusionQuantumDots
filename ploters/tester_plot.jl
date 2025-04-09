using CairoMakie
using DelimitedFiles

filename = "../build/evolution.dqmc.dat"
walkers_dist = readdlm(filename, comments = true)
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

x = collect(LinRange(xmin, xmax, nbins))

function psi_0(x)
    return @. (0.8 / π)^(0.25) * exp(-0.8 * x^2 / 2.)
end
##

with_theme(theme_latexfonts()) do
    fig = Figure();
    ax = Axis(fig[1,1], xlabel = "położenie [a.u.]", ylabel = "Ψ(x)",
              ylabelcolor = :blue, title = "Stan podstawowy oscylatora harmonicznego.");

    lines!(ax, x, walkers_dist[:, 1],
           label = "Rozkład wędrowców po czasie 1000 a.u.",
           color = :blue)

    lines!(ax, x, psi_0(x),
           label = "Rozwiązanie dokładne",
           color = :red,
           linestyle = :dash)

    display(fig)
end
