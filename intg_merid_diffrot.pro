;Paul A. Higgins - 10-mar-2012
;----------------------------------------------------------------------------->

function load_const, vcarr=vcarr, rsun=rsun, miss=miss

if keyword_set(vcarr) then retval=14.184 ;deg per sday
if keyword_set(rsun) then retval=6.955d8 ;in m
if keyword_set(miss) then retval=-1d6

return,retval
end

;----------------------------------------------------------------------------->
;TMM in decimal years
function load_hath, datfile, meridsmooth=msmth, drotsmooth=dsmth, yrminmax=yrmm

restore,datfile

miss=load_const(/miss)

;Do smoothing and interpolating beforehand to save processing time
fill_missing, MFARR, miss, 2
if n_elements(msmth) eq 2 then begin 
	MFARRsmth=smooth(MFARR,msmth,miss=miss,/edge)
	fill_missing, MFARRsmth, miss, 2
endif else MFARRsmth=MFARR

fill_missing, DRARR, miss, 2
if n_elements(dsmth) eq 2 then begin
	DRARRsmth=smooth(DRARR,dsmth,miss=miss,/edge)
	fill_missing, DRARRsmth, miss, 2
endif else DRARRsmth=DRARR

;Generate representative profiles for the time period of interest.
if n_elements(yrmm) eq 2 then begin
	wavg=where(DECYR ge yrmm[0] and DECYR le yrmm[1])
	mfarrtmm=average(MFARR[min(wavg):max(wavg),*],1,miss=miss)
	drarrtmm=average(DRARR[min(wavg):max(wavg),*],1,miss=miss)
endif else begin
	mfarrtmm=fltarr(n_elements(LATARR))
	drarrtmm=fltarr(n_elements(LATARR))
endelse

struct={LATARR:reform(LATARR), MFARR:MFARR, MFARRsmth:MFARRsmth, MFSDARR:MFSDARR, DRARR:DRARR, DRARRsmth:DRARRsmth, DRSDARR:DRSDARR, $
		CRS:CRS, AVGPROF:AVGPROF, AVGPROFIMG:AVGPROFIMG, AVGPROFMF:AVGPROFMF, AVGPROFMFIMG:AVGPROFMFIMG, DECYR:DECYR, $
		mfarrtmm:mfarrtmm, drarrtmm:drarrtmm}

return,struct
end

;----------------------------------------------------------------------------->
;convert a velocity profile from meters/s to degrees/s
function velmtodeg, LATARR, velarr, subcarr=subcarr

if keyword_set(subcarr) then begin
	;for differential rotation
	
	vdeg=(velarr)*360./(2.*!pi*load_const(/rsun)*cos(LATARR*!dtor))
	vdeg=vdeg+load_const(/vcarr)/3600./24.
	
endif else begin
	;for meridional flow (no cosine factor)
	vdeg=(velarr)*360./(2.*!pi*load_const(/rsun))
endelse

return,vdeg
end

;----------------------------------------------------------------------------->
;TIME=sec past 1-jan-1979 
;datfile=sav file with hathaway data
;deg=converts m/s prograde flow speed in ref frame of carrington coord 
;	to a rotation speed in degrees/s
pro load_merid, hathstr, time=thist, merid=outmerid, lat=outlat, tmatch=outtmatch, $
	boxsmooth=insmth, avg=avg, resid=resid, deg=deg, errstdv=outerr, tmm=tmm

;restore,datfile
;hathstr
MFARR=hathstr.MFARR
MFARRsmth=hathstr.MFARRsmth

miss=load_const(/miss)
if n_elements(insmth) eq 2 then MFARR=MFARRsmth ;MFARR=smooth(MFARR,insmth,miss=miss,/edge)
;fill_missing, MFARR, miss, 2


outlat=(hathstr.LATARR)[0,*]

;if TMM was set, uses the averaged profile from LOADHATH between the min and max times 
if keyword_set(tmm) then begin
	outmerid=hathstr.mfarrtmm
	if keyword_set(deg) then outmerid=velmtodeg(outlat, outmerid)
	if n_elements(insmth) eq 1 then outmerid=smooth(outmerid,insmth,miss=miss,/edge)
	return
endif

;Use profile that is averaged over solar cycle
if keyword_set(avg) then begin
	outmerid=hathstr.AVGPROFMF
	if keyword_set(deg) then outmerid=velmtodeg(outlat, outmerid)
	if n_elements(insmth) eq 1 then outmerid=smooth(outmerid,insmth,miss=miss,/edge)
	return
endif

thisdecyr=(thist)/3600./24./365.+1979.
wbest=where(abs(thisdecyr-hathstr.DECYR) eq min(abs(thisdecyr-hathstr.DECYR)))

outerr=(hathstr.MFSDARR)[wbest,*]
outtmatch=[((hathstr.DECYR)[wbest]-1979.)*3600.*24.*365.,(hathstr.DECYR)[wbest],(hathstr.CRS)[wbest]]

;Return the residual profile (this profile - SC averaged)
if keyword_set(resid) then begin
	outmerid=MFARR[wbest,*]-hathstr.AVGPROFMF
	if keyword_set(deg) then outmerid=velmtodeg(outlat, outmerid)
	if n_elements(insmth) eq 1 then outmerid=smooth(outmerid,insmth,miss=miss,/edge)
	return
endif

outmerid=(MFARR)[wbest,*]
if keyword_set(deg) then outmerid=velmtodeg(outlat, outmerid)

end

;----------------------------------------------------------------------------->

pro load_diffrot,hathstr, time=thist, diffrot=outdiffrot, lat=outlat, tmatch=outtmatch, $
	boxsmooth=insmth, avg=avg, resid=resid, deg=deg, errstdv=outerr, tmm=tmm

;restore,datfile
DRARR=hathstr.DRARR
DRARRsmth=hathstr.DRARRsmth

miss=load_const(/miss)
if n_elements(insmth) eq 2 then DRARR=DRARRsmth ;smooth(DRARR,insmth,miss=miss,/edge)
;fill_missing, DRARR, miss, 2
outlat=hathstr.LATARR[0,*]

if keyword_set(tmm) then begin
	outdiffrot=hathstr.drarrtmm
	if keyword_set(deg) then outdiffrot=velmtodeg(outlat, outdiffrot,/subcarr)
	if n_elements(insmth) eq 1 then outdiffrot=smooth(outdiffrot,insmth,miss=miss,/edge)
	return
endif

if keyword_set(avg) then begin
	outdiffrot=hathstr.AVGPROF
	if keyword_set(deg) then outdiffrot=velmtodeg(outlat, outdiffrot, /subcarr)
	if n_elements(insmth) eq 1 then outdiffrot=smooth(outdiffrot,insmth,miss=miss,/edge)
	return
endif

thisdecyr=(thist)/3600./24./365.+1979.
wbest=where(abs(thisdecyr-hathstr.DECYR) eq min(abs(thisdecyr-hathstr.DECYR)))

outerr=(hathstr.DRARR)[wbest,*]
outtmatch=[((hathstr.DECYR)[wbest]-1979.)*3600.*24.*365.,(hathstr.DECYR)[wbest],(hathstr.CRS)[wbest]]

if keyword_set(resid) then begin
	outdiffrot=DRARR[wbest,*]-hathstr.AVGPROF
	if keyword_set(deg) then outdiffrot=velmtodeg(outlat, outdiffrot, /subcarr)
	if n_elements(insmth) eq 1 then outdiffrot=smooth(outdiffrot,insmth,miss=miss,/edge)
	return
endif

outdiffrot=DRARR[wbest,*]
if keyword_set(deg) then outdiffrot=velmtodeg(outlat, outdiffrot, /subcarr)

end

;----------------------------------------------------------------------------->
;+
;NOTE: The Hathaway data goes to ~70. thus the max/min lat that will be used to iterate will be ~(+-70)
;Lats and lons in degrees
;ts in secs
;calculate the lon., lat., and t arrays for a point on the sun in HG coord. affected by the meridional flow and differential rotation
;Syntax: xx=intg_merid_diffrot(35.,0.,0.,flon=lon,flat=lat,ft=t,lonmax=0.)
;-
function intg_merid_diffrot, inilat, inilon, init, dt=indt, yrmm=inyrmm, $
	flon=outlon, flat=outlat, ft=outt, $
	lonmax=inlonmax, tmax=intmax, datfile=indatfile, $
	test=test, savprof=savprof, nomerid=nomerid, $
	hathstr=inhathstr, vmerid=invmerid, vdiffrot=invdiffrot

if n_elements(savprof) eq 1 then dosavprof=1 else dosavprof=0

if n_elements(indatfile) ne 1 then datfile='~/science/data/hathaway/hathaway_smart_compare.sav' else datfile=indatfile

if n_elements(inilat) ne 1 then ilat=40. else ilat=inilat
if n_elements(inilon) ne 1 then ilon=0. else ilon=inilon
if n_elements(init) ne 1 then it=0. else it=init

if n_elements(inyrmm) eq 0 then yrminmax=[1996.,1999.] else yrminmax=inyrmm

if n_elements(indt) ne 1 then dt=3600. else dt=indt ;s in one hr

;maximum longitude to reach and stop iterating (default=1 rotation)
if n_elements(inlonmax) ne 1 then lonmax=360. else lonmax=inlonmax
if n_elements(intmax) eq 1 then tmax=intmax

if n_elements(indt) ne 1 then dt=3600. else dt=indt

;smthmerid=[10,100.]
;smthdrot=[5,10]
smthmerid=[100.]
smthdrot=[10]

if n_elements(hathstr) eq 1 then hathstr=inhathstr else $
	hathstr=load_hath(datfile, meridsmooth=smthmerid, drotsmooth=smthdrot, yrminmax=[1996.,1999.])

thislat=ilat
thislon=ilon
thist=it
ltarr=ilat
lnarr=ilon
tarr=it

;!!!!use case statement to decide which variables to compare
;1. compare longitude with max longitude
;2. '' latitude '' max latitude
;3. '' time '' max time
;!!!

;TEMP!!??
if n_elements(invmerid) eq 0 or  n_elements(invdiffrot) eq 0 then begin
	load_merid, hathstr, time=thist, merid=vmerid, lat=mflatarr, boxsmooth=smthmerid, /deg, /tmm
	load_diffrot, hathstr, time=thist, diffrot=vdiffrot, lat=drlatarr, boxsmooth=smthdrot, /deg, /tmm
endif else begin
	vmerid=invmerid
	vdiffrot=invdiffrot
	mflatarr=reform((hathstr.LATARR)[0,*])
	drlatarr=reform((hathstr.LATARR)[0,*])
endelse

;Set iteration conditions
if n_elements(tmax) eq 1 then begin
	thisval=thist
	maxval=tmax
endif else begin
	thisval=thislon
	maxval=lonmax
endelse

;Iterate until the condition is met
while thisval lt maxval do begin

	;load the hathaway profiles for meridional flow and differential rotation
	;set tmm to use an average over some time range structure.mfarrtmm and drarrtmm
	;load_merid, hathstr, time=thist, merid=vmerid, lat=mflatarr, boxsmooth=smthmerid, /deg, /tmm
	;load_diffrot, hathstr, time=thist, diffrot=vdiffrot, lat=drlatarr, boxsmooth=smthdrot, /deg, /tmm

	;find best matching flow velocities to current position
	thismerid=vmerid[(where(abs(mflatarr-thislat) eq min(abs(mflatarr-thislat))))[0]]
	thisdiffrot=vdiffrot[(where(abs(drlatarr-thislat) eq min(abs(drlatarr-thislat))))[0]]
	
	if keyword_set(test) then print,thislat,thislon,thist

	;iterate to find next lat and lon and t
	if not keyword_set(nomerid) then $
		thislat=thislat+thismerid*dt
	thislon=thislon+thisdiffrot*dt
	thist=thist+dt

	;append new values to arrays
	ltarr=[ltarr,thislat]
	lnarr=[lnarr,thislon]
	tarr=[tarr,thist]

;Update iteration conditions
	if n_elements(tmax) eq 1 then thisval=thist $
		else thisval=thislon
endwhile

;SAVE THE EMPIRICAL PROFILES
latarr=hathstr.latarr
if dosavprof then save,vmerid,vdiffrot,latarr, file=savprof

outlat=ltarr
outlon=lnarr ;put in HG lon

 ;if value is past middle of backside of sun, shift to calc from 0 deg disk center going in solar east dir.

outt=tarr

thislon=(thislon mod 360.)
if thislon gt 180. then thislon=thislon-360.
retvar=[thislat,thislon,thist]
return, retvar

end

;----------------------------------------------------------------------------->