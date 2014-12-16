fdata = "sedov_taylor.txt"

n = read_four_arrays(fdata,r,rho,p,v)

plot,r,rho,xrange=[0,0.5],yrange=[0,4.2],xstyle=1,ystyle=1
end
