;Paul A. Higgins - 10-mar-2012
;----------------------------------------------------------------------------->
;Performs windowed smoothing on an array or image in the x-direction
;IMG		= Image to be smoothed
;WIDTH		= Width of smoothing window -> must be ODD!
;MISSVAL	= Value to be filled into regions that cannot be interpolated across
;				This value is also used to flag missing data when smoothing

function window_smooth, img, width=width, missval=missval, test=test

if n_elements(width) ne 1 then wdth=5. else wdth=width ;width of the window
if n_elements(missval) ne 1 then miss=-9999. else miss=missval ;missing value flag

dim=n_elements(size(img,/dim))
if wdth mod 2 eq 1 then wdth=wdth+1.

;determine window edges for each smoothed point
dwind=wdth-1. ;offset to use as pos rather than num.
x12arr=findgen(n_elements(img[*,0])-dwind)
x12arr=[transpose(x12arr),transpose(x12arr)]
x12arr[1,*]=x12arr[1,*]+dwind
x12arr=transpose(x12arr)

;run smoothing
imgsmth=img
for i=0.,(size(x12arr,/dim))[0]-1 do begin
	thiswind=imgsmth[x12arr[i,0]:x12arr[i,1],*]
	if dim eq 2 then avgslice=average(thiswind,1,miss=miss)
	if dim eq 1 then avgslice=average(thiswind,miss=miss)
	imgsmth[i+dwind/2.,*]=avgslice
endfor


if keyword_set(test) then begin
	if dim eq 1 then plot,imgsmth
	if dim eq 2 then plot_image,imgsmth
endif

return,imgsmth

end