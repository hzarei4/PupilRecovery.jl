using CoherentNoise, Chain
using ZernikePolynomials

using Plots, Printf
using View5D, FFTW
using PupilRecovery
using PSFDistiller

Base.show(io::IO, f::Float64) = @printf(io, "%.5f", f)

#sampler = opensimplex2_2d(seed=1)


# @time @chain opensimplex2_2d(seed=10) begin
#     fbm_fractal_2d(seed=10, source=_, frequency=5, persistence=0.7)
#     a1 = gen_image();
# end

# b1 = Float64.(Gray.(a1))

# t = mix(fbm_fractal_4d(seed=10, frequency=3, persistence=0.7), opensimplex2_4d(seed=10))
# y = gen_image(t)





function sampling_noise(
    sampler::S;
    w::Integer=400,
    h::Integer=400,
    xbounds::NTuple{2,Float64}=(-1.0, 1.0),
    ybounds::NTuple{2,Float64}=(-1.0, 1.0),
    time::Real
) where {N,S<:AbstractSampler{N}}
    x1, x2 = xbounds
    y1, y2 = ybounds
    xd = (x2 - x1) / w
    yd = (y2 - y1) / h
    img = Matrix{Float64}(undef, h, w)
    zw = ntuple(_->rand(sampler.random_state.rng, Float64) * 1000, N - 2) 

    Threads.@threads for x in 1:h
        cx = x * xd + x1
        for y in 1:w
            cy = y * yd + y1
            value = clamp(sample(sampler, cx, cy, 1.0, time) * 0.5 + 0.5, 0, 1)
            img[x, y] = value
        end
    end
    img
end


tmp = sampling_noise(
        fbm_fractal_4d(seed=10, frequency=2, persistence=0.1), w=100, h=100, time=1.0);
pupil_test = CreatePupil(phase = tmp, upsampling=4);

heatmap(angle.(pupil_test[25:75, 25:75]), 
        aspect_ratio=1.0, c=:twilight, clim=(-pi, pi), 
        legend = :none, axis=([], false)
        )

psf_test = CreatePSF(pupil_test)
@vt psf_test

psf_n_test = AddingNoise(psf_test, scaling=0.001) 
@vt psf_n_test


default(size=(800, 600))

@time x_solution, loss_trace = DMonPSF(psf_n_test; 
                        β=0.5, η=0.8, it_max=1000, tol=0.01, 
                        plotting=true, iterative_plotting=false,
                        );

heatmap(angle.(fftshift(x_solution)[50:150, 50:150]), 
        aspect_ratio=1.0, c=:twilight, clim=(-pi, pi), 
        legend = :none, axis=([], false)
        );


anim = @animate for i1 in 1:2

    begin
        pupil = CreatePupil(phase=(sampling_noise(fbm_fractal_4d(seed=10, frequency=4, persistence=0.2), time=i1/100.0)), upsampling=4)
        psf = CreatePSF(pupil)[150:250, 150:250]
        psf_n = AddingNoise(psf, scaling=0.1) 

        p00 = heatmap(angle.(pupil), 
                        title = "Simulated Phase", aspect_ratio=1.0, 
                        c=:twilight, clim=(-pi, pi), 
                        legend = :none,
                        );
        p01 = heatmap(psf_n, 
                        title = "Simulated PSF (cropped)", 
                        aspect_ratio=1, legend = :none,
                        #clim=(minimum(psf_n), maximum(psf_n)),
                        c= :grays,
                        );

        plot(p00, p01, layout=@layout([A B]), 
            framestyle=nothing, showaxis=false, 
            xticks=false, yticks=false, 
            size=(800, 400),  
            plot_title="Time=$(i1/50.0)",
            #plot_titlevspan=0.25
        )
    end


end;

gif(anim, fps=30)


gif(anim, "examples/noise_1.mp4", fps=60)