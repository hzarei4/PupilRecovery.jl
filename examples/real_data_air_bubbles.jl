
using Plots
using View5D, FFTW
using PupilRecovery
using PSFDistiller
using BioformatsLoader

BioformatsLoader.init()

beads_raw = bf_import("/home/hossein/Data/PSF/20220815_beads/AirBubbleStack_LowConc.czi", order="XYZCT")
beads=Float32.(beads_raw[1].data)[:,:,:,1,1]




main_img = beads[:,:,29]

@vt main_img



mypsf, rois, positions, selected, params, _ = distille_PSF(main_img, 
                                                roi_size=ntuple((x)->40,ndims(main_img)));
@vt mypsf

default(size=(900, 800))



x_solution, S_in = DMonPSF(mypsf, it_max=500, beta=0.2, plotting=true);


heatmap(abs.(abs.(fftshift(fft(x_solution))) .- abs.((mypsf)))./maximum(abs.(fftshift(mypsf))), 
title = "PSF diff", aspect_ratio=1)




@btime $x_solution = DMonPSF($mypsf, it_max=100, beta=0.5);

