
using Plots, Printf
using View5D, FFTW
using PupilRecovery
using PSFDistiller
using BioformatsLoader

Base.show(io::IO, f::Float64) = @printf(io, "%.8f", f)


BioformatsLoader.init()

beads_raw = bf_import("/home/hossein/Data/PSF/20220815_beads/AirBubbleStack_LowConc.czi", order="XYZCT")
beads=Float32.(beads_raw[1].data)[:,:,:,1,1]

@vv beads


main_img = beads[:,:,28]

@vt main_img



mypsf, rois, positions, selected, params, _ = distille_PSF(main_img, 
                                                roi_size=ntuple((x)->40,ndims(main_img)));
@vt mypsf

default(size=(1200, 800))



x_solution, loss_trace = DMonPSF(mypsf, it_max=1000, tol=0.01, beta=0.8, eta=0.8, plotting=true);



@btime $x_solution = DMonPSF($mypsf, it_max=100, beta=0.5);

