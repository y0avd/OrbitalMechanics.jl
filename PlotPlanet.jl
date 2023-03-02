X(r,theta,phi) = r * sin(theta) * sin(phi)
Y(r,theta,phi) = r * sin(theta) * cos(phi)
Z(r,theta,phi) = r * cos(theta)

function PlotPlanet(R,col)
    # spherical: (radius r, inclination θ, azimuth φ)
    thetas = range(0, stop=2pi,   length=50)
    phis   = range(-pi/2, stop=pi/2, length=50)

    xs = [X(R, theta, phi) for theta in thetas, phi in phis] 
    ys = [Y(R, theta, phi) for theta in thetas, phi in phis]
    zs = [Z(R, theta, phi) for theta in thetas, phi in phis]

    Plots.plot(xs, ys, zs, linetype=:surface, colorbar=false, color = col, fillalpha = 0.5);
end