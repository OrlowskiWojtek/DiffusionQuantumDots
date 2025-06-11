using CairoMakie
using DelimitedFiles

#filepath = "../build/results/ground_reblock"
#filepath = "
#outpath = "plots/ground_reblock.pdf"
#filepath = "../data/2el_1d_excited_10mev/reblock_analysis.dqmc.dat"
filepath = "../build/reblock_analysis.dqmc.dat"
data = readdlm(filepath, comments = true)[begin:end-1, :]

#

CairoMakie.activate!()

with_theme(theme_latexfonts()) do
    fig = Figure();
    ax_growth = Axis(fig[1,1], xlabel = "Iteracja reblokowania", ylabel = "Odchylenie standardowe σ estymatora wzrostu", ylabelcolor = :darkblue)

    ax_mixed = Axis(fig[1,1], yaxisposition = :right, ylabel = "Odchylenie standardowe σ estymatora mieszanego", ylabelcolor = :darkred)

    lines!(ax_mixed, data[:, 2], color = :darkred, label = "Błędy estymatora mieszanego")
    scatter!(ax_mixed, data[:, 2], color = :darkred)

    lines!(ax_growth, data[:, 3], color = :darkblue, label = "Błędy estymatora wzrostu")
    scatter!(ax_growth, data[:, 3], color = :darkblue)

    display(fig)
    #save("plots/1d_1el_excited_reblock.pdf", fig)
end
