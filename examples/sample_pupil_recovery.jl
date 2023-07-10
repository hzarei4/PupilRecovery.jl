using PupilRecovery
using Plots
using View5D
using FourierTools

upsampling = 1
pupil = CreatePupil(phase="zernike4", px=100, upsampling=upsampling)

# @vt angle.(pupil)


default(size=(400, 400))
heatmap(angle.(pupil), 
        aspect_ratio=1.0, c=:twilight, clim=(-pi, pi), 
        #legend = :none, 
        #axis=([], false),
        colorbar_title="Angle (Radians)",
        title="Random Phase with upsampling=$(upsampling)",
        )
# savefig("examples/phaseUpsample1.png")

psf = CreatePSF(pupil)

heatmap(psf, 
        title = "Simulated PSF with upsampling=$(upsampling)",
        aspect_ratio=1, 
        #legend = :none,
        colorbar_title=" ",
        #clim=(minimum(psf_n), maximum(psf_n)),
        c= :grays,
        )
# savefig("examples/psfUpsampling2.png")

psf_n = AddingNoise(psf, scaling=0.01)

heatmap(psf_n, 
        title = "Simulated PSF with upsampling=$(upsampling)
(Speckle pattern), Noise added", 
        aspect_ratio=1, 
        #legend = :none,
        colorbar_title=" ",
        #clim=(minimum(psf_n), maximum(psf_n)),
        c= :grays,
        )


heatmap(max.(0, sqrt.(psf_n) .- 1.5*mean(sqrt.(psf_n))), 
        title = "Simulated PSF with upsampling=$(upsampling)
(Speckle pattern), Noise added", 
        aspect_ratio=1, 
        #legend = :none,
        colorbar_title=" ",
        #clim=(minimum(psf_n), maximum(psf_n)),
        c= :grays,
        )
       
# @vt psf_n

# default(size=(900, 800))


# x_solution = DMonPSF(psf_n, it_max=200)
@time x_solution, loss_trace = DMonPSF(psf; β=0.99, η=0.8, it_max=2000, tol=1.0); #, plotting=true, iterative_plotting=false);

# savefig("examples/outputRandom.png")



begin
        
        p00 = heatmap(angle.(pupil), 
                title = "Simulated Phase", aspect_ratio=1.0, 
                c=:twilight, clim=(-pi, pi), 
                legend = :none,
                );
        p01 = heatmap(angle.(fftshift(x_solution)), 
                aspect_ratio=1.0, c=:twilight, clim=(-pi, pi), 
                legend = :none, 
                #axis=([], false),
                #colorbar_title="Angle (Radians)",
                title="Solution",
                );
        p02 = heatmap(angle.(fftshift(x_solution)) .- angle.(pupil), 
                aspect_ratio=1.0, c=:twilight, clim=(-pi, pi), 
                legend = :none, 
                #axis=([], false),
                #colorbar_title="Angle (Radians)",
                title="Difference",
                );
        plot(p00, p01, p02, layout=@layout([A B C]), 
        framestyle=nothing, showaxis=false, 
        xticks=false, yticks=false, 
        size=(1200, 300),  
        plot_title=" ",
        #plot_titlevspan=0.25
        )
end

savefig("examples/comparisonPhaseRandom.png")

