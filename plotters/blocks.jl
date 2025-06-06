using CairoMakie
using DelimitedFiles

#filepath = "../build/results/ground_reblock"
filepath = "../build/reblock_analysis.dqmc.dat"
#outpath = "plots/ground_reblock.pdf"

data = readdlm(filepath, comments = true)[begin:end-1, :]

# 

CairoMakie.activate!()

with_theme(theme_latexfonts()) do
    fig = Figure();
    ax = Axis(fig[1,1], xlabel = "Iteracja reblokowania", ylabel = "Odchylenie standardowe σ")

    lines!(ax, data[:, 2], color = :red, label = "Błędy estymatora mieszanego")
    scatter!(ax, data[:, 2], color = :red)

    lines!(ax, data[:, 3], color = :blue, label = "Błędy estymatora wzrostu")
    scatter!(ax, data[:, 3], color = :blue)

    Legend(fig[2,1], ax, framevisible = false, orientation = :horizontal)

    display(fig)
    save("plots/1d_1el_reblock.pdf", fig)
end
