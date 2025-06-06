using CairoMakie
using DelimitedFiles

filepath = "../build/results/ground_reblock"
outpath = "plots/ground_reblock.pdf"

data = readdlm(filepath, comments = true);

## 

CairoMakie.activate!()

with_theme(theme_latexfonts()) do
    fig = Figure();
    ax = Axis(fig[1,1], xlabel = "Iteracja reblokowania", ylabel = "Odchylenie standardowe σ")

    lines!(ax, data[:, 2], color = :blue)
    scatter!(ax, data[:, 2], color = :black)

    display(fig)
    #save(outpath, fig)
end
