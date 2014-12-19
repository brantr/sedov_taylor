fdata = "sedov_taylor.txt"

n = read_four_arrays(fdata,r,rho,p,v)

window,2,xsize=500,ysize=500
plot,r,rho,xrange=[0,0.5],yrange=[0,4.2],xstyle=1,ystyle=1


window,0,xsize=500,ysize=500
plot,r,p,xrange=[0,0.5],xstyle=1

end
