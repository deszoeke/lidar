#w_2pk_spectr.jl
#
# Simon's Julia script to show 2-peaked W-band velocity spectrum

# use the right python version
ENV["PYTHON"] = "/home/server/pi/homes/sdeszoek/anaconda/bin/python"

using DataArrays
using NetCDF
using PyPlot
# using DataFrames
# using Color
# using PyCall

filename="../../20083190000MMCRSpecMom.nc";
Heights=ncread(filename,"Heights",start=[1,1],count=[-1,1]);
time_offset=ncread(filename,"time_offset",start=[1],count=[-1]);
vNyq=ncread(filename,"NyquistVelocity",start=[1],count=[1]);
NumFFT=ncread(filename,"NumFFT",start=[1],count=[1]); # number of spectral estimates in a velocity spectrum
# lower bound of velocity bin, neg. vel. is upward away from radar
v=linspace(-vNyq[],vNyq[],NumFFT[]+1)
v=v[1:NumFFT[]]

spec=ncread(filename,"Spectra");
nheight=size(spec,2)
ntime=size(spec,3)

# set up colormap
cmap=get_cmap("RdYlBu_r")
bounds=[5.3:0.2:8.3] # discrete levels to color between
#norm=pycall(matplotlib[:colors][:BoundaryNorm],PyAny, bounds,cmap[:N]) # for PyCall v.1.1.1
norm=matplotlib[:colors][:BoundaryNorm](bounds, cmap[:N]) # for PyPlot v.1.4.3
# The syntax pythonlibrary[:sublib][:pythonfunction](arg1,arg2[:attrib],kwargs,etc)
# is the way of constructing Python objects (functions) in Julia.
# It often needs to be rewritten inside pycall:
# pycall(pythonlibrary[:sublib][:pythonfunction], PyAny, arg1,arg2[:attrib],kwargs,etc)

fit=findfirst(time_offset .>= 15*60)

clf()
ax=axes()
pcolormesh([v,vNyq],Heights[4:end]/1e3,log10(spec[:,4:end,fit]'),cmap=cmap,norm=norm); # uses too many colors
clim(bounds[[1,end]])
colorbar(spacing="uniform")
axis("image");
xylim=axis("tight")
axis([collect(xylim[1:2]), 0., 2.0])
ax[:set_xlim]([-vNyq[], vNyq[]]);
xlabel("Doppler velocity (m/s)")
ylabel("height (km)")
contour(v,Heights[4:end]/1e3,log10(spec[:,4:end,fit]'),[6.3, 6.9],colors="k")
savefig("spectr_eg_fig.eps")
savefig("spectr_eg_fig.png")

# spectral counts should never be negative, and they aren't (this hour)
if minimum(spec[:])<0
   println("Spectral counts less than 0 detected in $(filename).")
end

# test Julia's speed by doing big array operations for computing dynamic H-S noise floor
tic()
rtest=(1.0+NumFFT[])/NumFFT[]
noiseindex=zeros(Int16,nheight,ntime);
noisemean=zeros(Float32,nheight,ntime);
noisemax=zeros(Float32,nheight,ntime);

for itime=1:ntime
   specsort=sort(spec[:,:,itime],1)
   specsortsq=specsort.*specsort
   cumspec=cumsum(specsort)
   cumspecsq=cumsum(specsortsq)
   for iheight=1:nheight
      #noiseindex[iheight,itime]=maximum(find((cumspecsq[:,iheight].*[1:NumFFT[]]) .<= (rtest*cumspec[:,iheight].*cumspec[:,iheight])))
      noiseindex[iheight,itime]=findfirst((cumspecsq[:,iheight].*[1:NumFFT[]]) .> (rtest*cumspec[:,iheight].*cumspec[:,iheight]))-1
      noisemean[iheight,itime]=mean(specsort[1:noiseindex[iheight,itime],iheight],1)[]
      noisemax[iheight,itime]=specsort[noiseindex[iheight,itime],iheight]
   end
end
tdone=toc()

plot(Heights[4:end],noisemean[4:end,1:50:5000],".-")

allspecsort=sort(spec,1); # takes 12 s
pcolormesh(1:50:5000,Heights[4:end],log10(squeeze(allspecsort[65,4:end,1:50:5000],1)))

figure()
plot((60:70),spec[60:70,10:80,fit])


