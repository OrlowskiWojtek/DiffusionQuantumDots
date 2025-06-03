using CairoMakie
using DelimitedFiles

#filepath = "../build/results/ground_energies"
filepath = "../build/energies.dqmc.dat"
outpath = "plots/ground_energies.pdf"

data = readdlm(filepath, comments = true);
#

CairoMakie.activate!()

with_theme(theme_latexfonts()) do
    #outpath = "plots/ground_energies.pdf"
    steps = LinRange(1, length(data[:,1]), length(data[:,1]))

    fig = Figure();
    ax = Axis(fig[1,1], xlabel = "Krok czasowy", ylabel = "Średnia energia wędrowców")

    lines!(ax, data[:, 2], color = :red, label = "Estymator wzrostu")
    lines!(ax, data[:, 1], color = :black, label = "Estymator mieszany")
    #lines!(ax, steps, data, color = :black)
    #save(outpath, fig)
    axislegend()
    #outpath = "plots/ground_energies_with_blocks.pdf"
    #save(outpath, fig)
    display(fig)
end
