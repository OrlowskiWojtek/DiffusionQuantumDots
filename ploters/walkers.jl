using CairoMakie
using DelimitedFiles

filepath = "../build/results/ground_walkers"
outpath = "plots/ground_walkers.pdf"

walkers_dist = readdlm(filepath, comments = true)
ene = parse(Float64, split(readlines(filepath)[3], ":")[2])

function psi_0(x)
    return @. (0.4 / π)^(0.25) * exp(-0.4 * x^2 / 2.)
end

function V(x)
    return @. (0.16 * x^2)
end
##

with_theme(theme_latexfonts()) do
    fig = Figure();
    ax = Axis(fig[1:4,1:2], xlabel = "położenie [a.u.]", ylabel = "Ψ(x)",
              ylabelcolor = :blue, title = "Stan podstawowy oscylatora harmonicznego, energia = $ene");
    ax_r = Axis(fig[1:4,1:2], ylabel = "V(x)", yaxisposition = :right);
    hideydecorations!(ax_r)

    lines!(ax, walkers_dist[:,1], walkers_dist[:,2],
           label = "Rozkład wędrowców po czasie 1000 a.u.",
           color = :blue)

    lines!(ax, walkers_dist[:,1], psi_0(walkers_dist[:,1]),
           label = "Rozwiązanie dokładne",
           color = :red,
           linestyle = :dash)

    lines!(ax_r, walkers_dist[:,1], V(walkers_dist[:,1]),
           label = "Potencjał oscylatora harmonicznego",
           color = :black,
           linestyle = :dash)

    axes = [ax, ax_r]
    plots_in_fig = AbstractPlot[]
    labels_in_fig = AbstractString[]
    for ax in axes
        pl, lb = Makie.get_labeled_plots(ax, merge=false, unique=false)
        append!(plots_in_fig, pl)
        append!(labels_in_fig, lb)
    end

    ulabels = Base.unique(labels_in_fig)
    mergedplots = [[lp for (i, lp) in enumerate(plots_in_fig) if labels_in_fig[i] == ul]
        for ul in ulabels]

    Legend(fig[5,:], mergedplots, ulabels, framevisible = false)

    save(outpath, fig)
    display(fig)
end
