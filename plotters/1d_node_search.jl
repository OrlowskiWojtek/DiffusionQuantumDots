using CairoMakie
using DelimitedFiles
using LsqFit

data = readdlm("../data/1el_1d/excited_basic/zero_place_search", comments = true)

model(x, p) =@. x * p[1] + p[2]
left_data = data[begin:8,:]
right_data = data[8:end,:]
p0 = [0.5, 0.5]

left_fit = curve_fit(model, left_data[:,1], left_data[:,2], p0)
right_fit = curve_fit(model, right_data[:,1], right_data[:,2], p0)

left_line = model(left_data[:,1], left_fit.param)
right_line = model(right_data[:,1], right_fit.param)
##

set_theme!(theme_latexfonts())

fig = Figure();
ax = Axis(fig[1,1], xlabel = "Położenie węzła [nm]", ylabel = "Estymator wzrostu [meV]");

lines!(ax, right_data[:,1] .* 0.059, right_line, color = :black, label = "Dopasowane proste")
lines!(ax, left_data[:,1] .* 0.059, left_line, color = :black, label = "Dopasowane proste")
scatter!(ax, data[:,1] * 0.059, data[:,2], color = :darkred, label = "Estymator wzrostu")
errorbars!(ax, data[:,1] * 0.059, data[:,2], data[:,3], color = :darkred, whiskerwidth = 10)
#Legend(fig[2,1], orientation = :horizontal, unique = true)
axislegend(unique = true)

save("plots/1el_1d_node_error.pdf", fig)
display(fig);

##
left_fit.param[1] / 0.059
right_fit.param[1] / 0.059

