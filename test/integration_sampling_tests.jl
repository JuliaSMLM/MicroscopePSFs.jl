let
    function reference_integral(func, camera, emitter, sampling)
        nx = length(camera.pixel_edges_x) - 1
        ny = length(camera.pixel_edges_y) - 1
        result = zeros(typeof(emitter.photons), ny, nx)

        for ix in 1:nx, iy in 1:ny
            x_start = camera.pixel_edges_x[ix]
            x_end = camera.pixel_edges_x[ix + 1]
            y_start = camera.pixel_edges_y[iy]
            y_end = camera.pixel_edges_y[iy + 1]
            dx = x_end - x_start
            dy = y_end - y_start

            if sampling == 1
                x = (x_start + x_end)/2
                y = (y_start + y_end)/2
                value = func(x - emitter.x, y - emitter.y)
                result[iy, ix] = value * (dx * dy)
            else
                dx_sample = dx / sampling
                dy_sample = dy / sampling
                area_factor = dx_sample * dy_sample
                xs = range(x_start + dx_sample/2, x_end - dx_sample/2, length=sampling)
                ys = range(y_start + dy_sample/2, y_end - dy_sample/2, length=sampling)
                first_val = func(xs[1] - emitter.x, ys[1] - emitter.y)
                pixel_sum = zero(typeof(first_val))
                for x in xs, y in ys
                    pixel_sum += func(x - emitter.x, y - emitter.y)
                end
                result[iy, ix] = pixel_sum * area_factor
            end
        end

        result .*= emitter.photons
        return result
    end

    cameras = (
        IdealCamera(collect(0.0:0.11:0.44), collect(0.0:0.09:0.27)),
        IdealCamera([-0.03, 0.07, 0.19, 0.34, 0.5], [-0.02, 0.09, 0.21, 0.37]),
    )
    airy = AiryPSF(1.2, 0.68)
    scalar = ScalarPSF(1.2, 0.68, 1.33; pupil_data=ones(ComplexF64, 32, 32))
    vector = VectorPSF(1.2, 0.68, DipoleVector(0.6, 0.0, 0.8);
        n_medium=1.33, n_coverslip=1.33, n_immersion=1.33, grid_size=32)
    fixtures = (
        (airy, Emitter2D(0.18, 0.14, 1.0)),
        (scalar, Emitter3D(0.18, 0.14, 0.12, 1.0)),
        (vector, Emitter3D(0.18, 0.14, 0.12, 1.0)),
    )

    @testset "Independent range reference" begin
        for (psf, emitter) in fixtures
            func = if emitter isa Emitter2D
                (x, y) -> psf(x, y)
            else
                (x, y) -> psf(x, y, emitter.z)
            end
            for camera in cameras, sampling in (1, 2, 3, 4)
                expected = reference_integral(func, camera, emitter, sampling)
                for threaded in (false, true)
                    actual = integrate_pixels(psf, camera, emitter;
                        sampling=sampling, threaded=threaded)
                    @test isequal(actual, expected)
                end
            end
        end
    end

    @testset "Two-point endpoints" begin
        for T in (Float32, Float64), camera in cameras
            for edges in (T.(camera.pixel_edges_x), T.(camera.pixel_edges_y))
                for i in 1:(length(edges) - 1)
                    first_edge, last_edge = edges[i], edges[i + 1]
                    step = (last_edge - first_edge) / 2
                    first_point = first_edge + step/2
                    last_point = last_edge - step/2
                    points = range(first_point, last_point, length=2)
                    @test Tuple(points) === (first_point, last_point)
                end
            end
        end
    end
end
