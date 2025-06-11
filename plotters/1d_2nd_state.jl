using CairoMakie
using DelimitedFiles
using LsqFit

data = readdlm("../data/1el_1d/excited_basic/second_excited_search", comments = true)

##

set_theme!(theme_latexfonts())

fig = Figure();
ax = Axis(fig[1,1], xlabel = "Położenie węzła [nm]", ylabel = "Estymator wzrostu [meV]");

scatter!(ax, data[:,1] * 0.059, data[:,2], color = :darkred, label = "Estymator wzrostu")
errorbars!(ax, data[:,1] * 0.059, data[:,2], data[:,3], color = :darkred, whiskerwidth = 10)
axislegend(unique = true)

save("plots/1el_1d_nodes_2nd_excited.pdf", fig)
display(fig);

