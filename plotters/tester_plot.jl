using CairoMakie
using DelimitedFiles

CairoMakie.activate!()

a = 20.
b = 5.
jastrow(r) = @. exp(-1. / 2. * a / r * (1. - exp(-r / b)))

r_vals = LinRange(0.1, 5 / 0.0529, 50)
vals = jastrow(r_vals)

fig = Figure();
ax = Axis(fig[1,1]);

lines!(ax, r_vals, vals);

display(fig)
