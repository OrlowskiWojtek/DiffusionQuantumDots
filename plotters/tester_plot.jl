using CairoMakie
using DelimitedFiles

CairoMakie.activate!()

function ground_state_2d_harmonic_oscillator(x, y; homegax= 3.0 / 27211.6, homegay= 5.0 / 27211.6)
    # Funkcja falowa stanu podstawowego dla oscylatora 1D w wymiarze x
    psi_x(x) = (homegax / π)^(1/4) * exp(-0.5 * homegax * x^2)
    
    # Funkcja falowa stanu podstawowego dla oscylatora 1D w wymiarze y
    psi_y(y) = (homegay / π)^(1/4) * exp(-0.5 * homegay * y^2)
    
    # Stan podstawowy w 2D: iloczyn stanów podstawowych x i y
    return psi_x(x) * psi_y(y)
end
x = [_x for _x in -20:0.5:20] ./ 0.0529
y = [_y for _y in -20:0.5:20] ./ 0.0529

psi = [ground_state_2d_harmonic_oscillator(_x, _y) for _x in x, _y in y]

fig = Figure();
ax = Axis3(fig[1,1], xlabel = "x", ylabel = "y")

surface!(ax, x, y, psi)

display(fig)
