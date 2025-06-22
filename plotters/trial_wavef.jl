using GLMakie, CairoMakie
using DelimitedFiles

filename = "../data/2el_1d_excited_10mev_guided/TrialWavefunctionTest"
data = readdlm(filename, comments = true)

file = readlines(filename)
xmin = -50
xmax = 50 
nbins = size(data, 1)

x = y = LinRange(xmin, xmax, nbins)
#

CairoMakie.activate!()

with_theme(theme_latexfonts()) do
    fig = Figure(dpi = 300);
    ax = Axis3(fig[1,1], 
               xlabel = "x1 [nm]",
               ylabel = "x2 [nm]")

    hidezdecorations!(ax)

    cm = surface!(ax, x, y, data,
             colormap = :coolwarm,
             label = "Ψ(x,y)")

    ticks_psi = [ round(val, digits = 3) for val in collect(LinRange(findmin(data)[1], findmax(data)[1], 6)) ]

    Colorbar(fig[2,1], cm, label = L"$\Psi$(\mathbf{R})", vertical = false)

    display(fig)
    save("plots/10mev_excited_trial_wavef.png", fig)
end

## check for correct normalisation <- it is ok

dx = x[2] - x[1]
s = sum(data.^2) * dx^2


