pro ferengi_analysis


;--------------------------------------------------------------------
; ## Set up plotting
  cleanplot, /silent
  set_plot, 'ps'
  device, filename='fake_results.ps', $
          /color, /landscape, bits_per_pixel=8,  $
          /helvetica
  loadct, 0, /silent
  set_thick,4

  !p.charsize=2
  !p.font=0
;--------------------------------------------------------------------



; Note: We are testing on: T01_SMOOTH_OR_FEATURES_A02_FEATURES_FRAC
; for now.





;--------------------------------------------------------------------
; # Read in data and required metadata
; Read in Brooke's first most certainly incorrect reduction FITS file
data=mrdfits('../../data/ferengi_classifications_collated.fits',1)

; Read in Edmond's metadata file
info=mrdfits('../../data/GZ2_FERENGI_objects_FINAL_7_6_13.fits',1)

; Read in Stuart's metadata file, converted to FITS via TOPCAT
meta=mrdfits('../../data/galaxy_zoo_ferengi_subjects.fits',1)
;--------------------------------------------------------------------





  plot, [0], [0], $
        xr=[0, 1.2], xstyle=1, xtitle='redshift', $
        yr=[0, 1], ystyle=1, ytitle=textoidl('vote fraction for question X, f_{x,z}'), $
        xtickv=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2], xticks=6, $
        ytickv=[0.0, 0.25, 0.50, 0.75, 1.0], yticks=5


; Redshift array
  redshift_array_1=[0.030, 0.3+findgen(8)*0.1]


; Loop over each galaxy in Edmond's catalog
for i=0L,n_elements(info.objid)-1 do begin
                                ; Print ID and some info
   print, info[i].objid
   
                                ; Match the current object in
                                ; Stuart's table. 
   ii=where(meta.sdss_id eq strtrim(string(info[i].objid),2))


   stop
endfor


;--------------------------------------------------------------------
; ## Clean up
  device, /close
  set_plot, 'x'
  cleanplot, /silent
;--------------------------------------------------------------------

end
