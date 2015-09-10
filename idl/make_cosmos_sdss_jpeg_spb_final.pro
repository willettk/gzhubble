;+
; NAME:
;   make_hdf_jpeg
; BUGS:
;   everything hard wired
;-
pro make_cosmos_sdss_jpeg_spb_final,scale=scale,nonlin=nonlin,objno=objno

;read in catalog

;y=mrdfits('/Users/rgriffit/ACS/catalogs/cosmos_i_public_catalog_4.0.fits',1)

;b=where(y.mag_best_hi le 24.5 and y.mag_best_hi ge 10 and y.objno eq objno)
;b=where(y.mag_best_hi le 24.5 and y.mag_best_hi ge 10)

;y=y[b]

;for i=0l,n_elements(y)-1 do begin

for i=0l,1 do begin
;imgr='/Volumes/ACS/COSMOS/djsin/'+strcompress(y[i].objno,/rem)+'.r.fits'
;imgg='/Volumes/ACS/COSMOS/djsin/'+strcompress(y[i].objno,/rem)+'.g.fits'
;imgb='/Volumes/ACS/COSMOS/djsin/'+strcompress(y[i].objno,/rem)+'.b.fits'

;imgr=''+strcompress(y[i].objno,/rem)+'.r.fits'
;imgg=''+strcompress(y[i].objno,/rem)+'.g.fits'
;imgb=''+strcompress(y[i].objno,/rem)+'.b.fits'

imgr=''+strcompress(objno,/rem)+'.r.fits'
imgg=''+strcompress(objno,/rem)+'.g.fits'
imgb=''+strcompress(objno,/rem)+'.b.fits'

;output=''+strcompress(y[i].objno,/rem)+'_'+strcompress(scale,/rem)+'_'+strcompress(nonlin,/rem)+''
output=''+strcompress(objno,/rem)+'_'+strcompress(scale,/rem)+'_'+strcompress(nonlin,/rem)+''
;output='/Users/rgriffit/Desktop/COSMOS_zoo_images/'+strcompress(y[i].objno,/rem)+''

; read data
imgr=mrdfits(imgr,0,hdi)
imgb=mrdfits(imgb,0,hdv)
imgg=(imgr+imgb)/2

fscale = scale * 1.0d-4 * [0.8, 1.0, 1.0]

djs_rgb_make,imgr,imgg,imgb,scale=fscale,nonlinearity=nonlin,satvalue=0,name=output+'.tiff',/tiff

endfor

return

end
