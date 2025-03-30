using GLMakie, CairoMakie
using DelimitedFiles

filename = "../build/results/1st_state_2d_anharmonic"
data = readdlm(filename, comments = true)

x = y = LinRange(-5., 5., 200)

dx = x[2] - x[1]

##

GLMakie.activate!()

with_theme(theme_latexfonts()) do
    fig = Figure();
    ax = Axis3(fig[1,1], 
               xlabel = "x [a.u]",
               ylabel = "y [a.u]",
               zlabel = "Ψ(x,y)")

    surface!(ax, x, y, data,
             colormap = :bluesreds)

    display(fig)
end

## sanity normalization check

norm = sum(data.^2) * dx^2

