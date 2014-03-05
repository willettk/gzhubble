pro ferengi_analysis


;--------------------------------------------------------------------
; ## Set up plotting
  cleanplot, /silent
  set_plot, 'ps'
  device, filename='somewhat_fake_results.ps', $
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

; Read in the GZ2 classifications for the original galaxies, generated
;by Brooke
gz2=mrdfits('../../data/GZ2_FERENGI_matched_to_gz2_catalog.fits',1)
;--------------------------------------------------------------------




; Redshift array
  redshift_array_1=[0.030, 0.3+findgen(8)*0.1]

; Array for storing the vote fraction and redshift of the original
; galaxies
original_galaxy_vote     = dblarr(n_elements(info.objid))
original_galaxy_redshift = dblarr(n_elements(info.objid))


; Loop over each galaxy in Edmond's catalog
for i=0L,n_elements(info.objid)-1 do begin
                                ; Print ID and some info
   print, info[i].objid
   ;print, info[i].desired_z
   ;print, info[i].sb_bin
   
                                ; Process the original GZ2
                                ; classification of the current
                                ; object. No matching required as
                                ; Brooke simply matched
                                ; Edmond's table (info) to the
                                ; GZ2 database

   ;print, gz2[i].T01_SMOOTH_OR_FEATURES_A02_FEATURES_OR_DISK_WEIGHTED_FRACTION
   original_galaxy_vote[i]    = gz2[i].T01_SMOOTH_OR_FEATURES_A02_FEATURES_OR_DISK_WEIGHTED_FRACTION
   

                                ; Match the current object in
                                ; Stuart's table. 
   ii=where(meta.sdss_id eq strtrim(string(info[i].objid),2))
   original_galaxy_redshift[i] = meta[ii[0]].redshift
   ;print, original_galaxy_redshift

                                ; ### Experimental part. I'm
                                ; trying to reverse-engineer
                                ; Edmond's scheme to locate
                                ; results for the FERENGI'd galaxies
   
                                ; Loop over ALL the
                                ; subject_id's that correspond
                                ; to the current galaxy
                                ; according to Stuart's table
   n_ferengi=n_elements(ii)
   for j=0L,n_elements(ii)-1 do begin
      print, meta[ii[j]].subject_id
      current_subject_id = meta[ii[j]].subject_id

      ; Locate the current subject in Brooke's data
      jj=where(data.subject_id eq current_subject_id)
      
      ; Extract the vote fraction
      current_subject_id_vote = data[jj].T01_SMOOTH_OR_FEATURES_A02_FEATURES_FRAC
      ;print, current_subject_id_vote

      
   endfor
   stop

endfor


  plot, [0], [0], $
        xr=[0, 1.1], xstyle=1, xtitle='redshift', $
        yr=[0, 1], ystyle=1, ytitle=textoidl('vote fraction for question X, f_{x,z}'), $
        xtickv=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2], xticks=6, $
        ytickv=[0.0, 0.25, 0.50, 0.75, 1.0], yticks=5

  oplot, original_galaxy_redshift, original_galaxy_vote, psym=sym(1), symsize=0.8


stop

;--------------------------------------------------------------------
; ## Clean up
  device, /close
  set_plot, 'x'
  cleanplot, /silent
;--------------------------------------------------------------------

end
