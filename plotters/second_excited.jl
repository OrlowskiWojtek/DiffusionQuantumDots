using CairoMakie
using DelimitedFiles

#filename = "../data/1el_1d/excited_basic/evolution.dqmc.dat"
filename = "../build/evolution.dqmc.dat"
data = readdlm(filename, comments = true)

psi = data[:, 1]

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

x = LinRange(xmin, xmax, nbins)

x_atomic = x ./ 0.0529

function psi_2(x, homega, mass)
    α = homega * mass
    return @. (1. / sqrt(2^2 * factorial(2)) * (α / π)^(0.25) * (4 * α * x^2 - 2) * exp(-α * x^2 / 2.))
end

psi_exact = abs.(psi_2(x_atomic, 3 / 27211.6, 0.067))

#

with_theme(theme_latexfonts()) do
    fig = Figure();
    ax = Axis(fig[1,1], 
              xlabel = "x [nm]",
              ylabel = "|Ψ(x)|")

    lines!(ax, x, psi_exact, color = :red, label = "Rozwiązanie dokładne")
    scatter!(ax, x, psi, color = :black, label = "Rozwiązanie DMC")
 
    axislegend()
    display(fig)
    save("plots/1d_1el_2nd_excited.pdf", fig)
end
