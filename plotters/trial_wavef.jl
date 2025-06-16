using GLMakie, CairoMakie
using DelimitedFiles

filename = "../build/TrialWavefunctionTest"
data = readdlm(filename, comments = true)

file = readlines(filename)
xmin = -50
xmax = 50 
nbins = size(data, 1)

x = y = LinRange(xmin, xmax, nbins)
#

GLMakie.activate!()

with_theme(theme_latexfonts()) do
    fig = Figure();
    ax = Axis3(fig[1,1], 
               xlabel = "x [nm]",
               ylabel = "y [nm]")

    hidezdecorations!(ax)

    cm = surface!(ax, x, y, abs.(data),
             colormap = :coolwarm,
             label = "Ψ(x,y)")

    ticks_psi = [ round(val, digits = 3) for val in collect(LinRange(findmin(data)[1], findmax(data)[1], 6)) ]

    Colorbar(fig[2,1], cm, label = "Ψ(x,y)", vertical = false)

    #xlims!(ax, (xmin, xmax))
    #ylims!(ax, (xmin, xmax))
    zlims!(ax, (0, nothing))

    display(fig)
    #save("plots/example_trial_wavef.pdf", fig)
end

## check for correct normalisation <- it is ok

dx = x[2] - x[1]
s = sum(data.^2) * dx^2


