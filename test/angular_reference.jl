using LinearAlgebra

# Independent angular integration for matched indices. No pupil/field kernels are used.
# Fresnel transmission is 1; apodization is 1/sqrt(cos(theta)). The area element
# is (n/lambda)^2*sin(theta)*cos(theta) dtheta dphi. Energy is integrated analytically.
function part5_angular_reference(NA, n, lambda, points; order=128, azimuths=256)
    beta = [i/sqrt(4*i^2-1) for i in 1:order-1]
    rule = eigen(SymTridiagonal(zeros(order), beta))
    theta_max = asin(NA/n)
    theta = (rule.values .+ 1) .* (theta_max/2)
    weights = theta_max .* rule.vectors[1,:].^2
    amplitudes = [zeros(ComplexF64,2,3) for _ in points]
    k = n/lambda
    for (t,w) in zip(theta,weights)
        s,c = sincos(t)
        measure = w*(2pi/azimuths)*k^2*s*sqrt(c)
        for j in 0:azimuths-1
            sp,cp = sincos(2pi*j/azimuths)
            ex = (c*cp^2+sp^2, (c-1)*sp*cp, s*cp)
            ey = ((c-1)*sp*cp, c*sp^2+cp^2, s*sp)
            for (point,a) in zip(points,amplitudes)
                x,y,z = point
                ph = measure*cis(2pi*k*(s*cp*x+s*sp*y+c*z))
                for d in 1:3
                    a[1,d] += ph*ex[d]
                    a[2,d] += ph*ey[d]
                end
            end
        end
    end
    c = sqrt(1-(NA/n)^2)
    transverse_energy = pi*k^2*((1-c)+(1-c^3)/3)
    axial_energy = 2pi*k^2*((1-c)-(1-c^3)/3)
    energy = 2transverse_energy + axial_energy
    d = [0.6,0.0,0.8]
    fixed_energy = 0.6^2*transverse_energy + 0.8^2*axial_energy
    rotating = [sum(abs2,a)/energy for a in amplitudes]
    fixed = [sum(abs2,a*d)/fixed_energy for a in amplitudes]
    return rotating, fixed
end
