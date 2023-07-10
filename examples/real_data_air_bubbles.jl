
using Plots, Printf
using View5D, FFTW
using PupilRecovery
using PSFDistiller
using BioformatsLoader

Base.show(io::IO, f::Float64) = @printf(io, "%.5f", f)


BioformatsLoader.init()

beads_raw = bf_import("/home/hossein/Data/PSF/20220815_beads/AirBubbleStack_LowConc.czi", order="XYZCT")
beads=Float32.(beads_raw[1].data)[:,:,:,1,1]

@vv beads


main_img = beads[:,:,26]

@vt main_img



mypsf, rois, positions, selected, params, _ = distille_PSF(main_img, 
                                                roi_size=ntuple((x)->40,ndims(main_img)));
@vt mypsf

default(size=(900, 600))


@time x_solution, loss_trace = DMonPSF(mypsf; β=0.99, η=0.7, it_max=1000, tol=0.001, plotting=false, iterative_plotting=false);



savefig("examples/avgBeadsResult.png")


rois_all_phase = zeros(size(mypsf)..., length(rois))

for i in 1:length(rois)

    @time x_solution, loss_trace = DMonPSF(rois[i]; β=0.99, η=0.7, it_max=1000, tol=0.015, plotting=false, iterative_plotting=false);
    rois_all_phase[:, :, i] = angle.(x_solution)
    #sleep(2)
end

heatmap((fftshift(rois_all_phase[:, :, 2])), 
        aspect_ratio=1.0, c=:twilight, clim=(-pi, pi), 
        legend = :none, axis=([], false)
        )


anim = @animate for i1 in 1:length(rois)
    begin

        p00 = heatmap(fftshift(rois_all_phase[:, :, i1]), 
                        title = "Recovered Phase", aspect_ratio=1.0, 
                        c=:twilight, clim=(-pi, pi), 
                        legend = :none,
                        );
        p01 = heatmap(rois[i1], 
                        title = "Measured PSF", 
                        aspect_ratio=1, legend = :none,
                        #clim=(minimum(psf_n), maximum(psf_n)),
                        c= :grays,
                        );
        
        p02 = heatmap(log10.(main_img), 
                        title = "Main Image (focus point)", 
                        aspect_ratio=1, legend = :none,
                        #clim=(minimum(psf_n), maximum(psf_n)),
                        c= :grayC,
                        );
        quiver!([positions[i1][2] + 50], [positions[i1][1] + 30], quiver=([-45.0], [-25.0]), c="black", lw=5)
        annotate!([positions[i1][2] + 70], [positions[i1][1] + 50],"$")

        plot(p00, p01, p02, layout=@layout([A B C]), 
            framestyle=nothing, showaxis=false, 
            xticks=false, yticks=false, 
            size=(2000, 700),  
            plot_title="$(i1)",
            #plot_titlevspan=0.25
        )

    end


end;
    
gif(anim, fps=0.5)
        
        
gif(anim, "examples/rois_recovered.mp4", fps=0.5)