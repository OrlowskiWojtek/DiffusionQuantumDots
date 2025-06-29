using CairoMakie
using DelimitedFiles

filename = "../build/evolution189.035917.dqmc.dat"
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

function psi_0(x, homega, mass)
    return @. (homega * mass/ π)^(0.25) * exp(-homega * mass * x^2 / 2.)
end

psi_exact = psi_0(x_atomic, 3 / 27211.6, 0.067)
##
#

with_theme(theme_latexfonts()) do
    fig = Figure();
    ax = Axis(fig[1,1], 
              xlabel = "x [nm]",
              ylabel = "Ψ(x)")

    #lines!(ax, x, psi_exact, color = :red, label = "Rozwiązanie dokładne")
    scatter!(ax, x, psi, color = :black, label = "Rozwiązanie DMC")
 
    axislegend()
    display(fig)
    #save("plots/1d_1el_ground.pdf", fig)
end

##

dx = x[2] - x[1]
s = sum(psi.^2) * (dx / 0.0529)

print(sum(psi_exact.^2) * (dx / 0.0529))

