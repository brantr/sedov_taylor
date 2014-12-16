function read_four_arrays,fdata,x,y,z,w

	openr,1,fdata
	n = 0L
	readf,1,n
	x = dblarr(n)
	y = dblarr(n)
	z = dblarr(n)
	w = dblarr(n)
	xb = 0.0D
	yb = 0.0D
	zb = 0.0D
	wb = 0.0D
	for i=0L,n-1 do begin
		readf,1,xb,yb,zb,wb
		x(i) = double(xb)
		y(i) = double(yb)
		z(i) = double(zb)
		w(i) = double(wb)
	endfor
	close,1
	return,n


end
