;predict the disk-crossings for given seed points
;use to select mosaic images
pro predict_dcross


;load in the two seed points...


;[time,lat,lon] vs seed index
nseed=2
seeds=fltarr(nseed,1,3)
seeds[*,0,0]=[5.9054400d8,5.9140800d8] ;time
seeds[*,0,1]=[30.4416,33.9541] ;lat
seeds[*,0,2]=[0.,0.] ;lon

ncross=6

;array of seed index vs crossing index vs [time, lat, lon]
keypts=fltarr(nseed,ncross,3)
keypts[*,0,*]=seeds



for j=0,nseed-1 do begin


	for i=1,ncross-1 do begin
		lastkey=keypts[j,i-1.,*]
		;print,lastkey
	
		outparam=intg_merid_diffrot(lastkey[1], lastkey[2], lastkey[0], dt=3600., $
			flon=outlon, flat=outlat, ft=outt, $
			lonmax=360.) ;,datfile=indatfile, test=test, savprof=savprof, nomerid=nomerid
		
		print,outparam
		
		;[thislat,thislon,thist]
		keypts[j,i,*]=[outparam[2],outparam[0],outparam[1]]		
		
		
		
		
	endfor


endfor

print,anytim(keypts[0,*,0],/vms)

;first set:
set0=['18-Sep-1997 00:00:00.000',$
	'14-Oct-1997 01:12:32.000',$
	'9-Nov-1997 06:24:00.000',$
	'5-Dec-1997 15:34:24.000',$
	'1-Jan-1998 06:43:12.000',$
	'28-Jan-1998 03:50:24.000']

;second set:
set1=['28-Sep-1997 00:00:00.000',$
	'24-Oct-1997 09:10:24.000',$
	'20-Nov-1997 00:19:12.000',$
	'16-Dec-1997 20:26:40.000',$
	'12-Jan-1998 23:32:16.000',$
	'9-Feb-1998 08:36:16.000']

save,set0,set1,file=path+'mosaic_times.sav',/ver
stop

end