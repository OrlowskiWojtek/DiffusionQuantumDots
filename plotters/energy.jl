using CairoMakie
using DelimitedFiles

filepath = "../build/results/ground_energies"
outpath = "plots/ground_energies.pdf"

data = readdlm(filepath, comments = true);

## 
with_theme(theme_latexfonts()) do
    outpath = "plots/ground_energies.pdf"
    steps = LinRange(1, length(data[:,1]), length(data[:,1]))

    fig = Figure();
    ax = Axis(fig[1,1], xlabel = "Krok czasowy", ylabel = "Średnia energia wędrowców")

    #lines!(ax, data[begin:1000:end, 2], color = :black)
    lines!(ax, steps[begin:1000:end], data[begin:1000:end, 1], color = :black)

    block_size = 2^16
    blocks = [block_size * i for i in 1:((length(steps)) / block_size)]

    save(outpath, fig)
    vlines!(ax, blocks, color = :blue, label = "Podział na bloki")
    xlims!(ax, (1e6, 3e6))

    axislegend()
    outpath = "plots/ground_energies_with_blocks.pdf"
    save(outpath, fig)
    display(fig)
end
