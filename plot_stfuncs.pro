fdata = "stfuncs.txt"

n = read_four_arrays(fdata,xi,G,Z,V)

window,0,xsize=400,ysize=400
plot,xi,G

window,1,xsize=400,ysize=400
plot,xi,Z

window,2,xsize=400,ysize=400
plot,xi,V
end
