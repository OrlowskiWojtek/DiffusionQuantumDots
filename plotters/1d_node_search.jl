using CairoMakie
using DelimitedFiles

data = readdlm("../data/1el_1d/excited_basic/zero_place_search", comments = true)

##

fig = Figure();
ax = Axis(fig[1,1], xlabel = "Położenie węzła [nm]", ylabel = "Estymator wzrostu [meV]");

lines!(ax, data[:,1] * 0.059, data[:,2], color = :black)
scatter!(ax, data[:,1] * 0.059, data[:,2], color = :red)

display(fig);
