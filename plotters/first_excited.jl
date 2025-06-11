using CairoMakie
using DelimitedFiles

#filename = "../build/evolution.dqmc.dat"
filename = "../data/1el_1d/excited_basic/evolution.dqmc.dat"
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

function psi_1(x, homega, mass)
    return @. abs(1. / sqrt(2) * (homega * mass / π)^(0.25) * exp(-homega * mass * x^2 / 2.) * sqrt(homega * mass) * x * 2)
end

psi_exact = psi_1(x_atomic, 3 / 27211.6, 0.067)

##

with_theme(theme_latexfonts()) do
    fig = Figure();
    ax = Axis(fig[1,1], 
              xlabel = "x [nm]",
              ylabel = "|Ψ(x)|")

    lines!(ax, x, psi_exact, color = :red, label = "Rozwiązanie dokładne")
    scatter!(ax, x, psi, color = :black, label = "Rozwiązanie DMC")
 
    axislegend()
    display(fig)
    save("plots/1d_1el_excited.pdf", fig)
end
