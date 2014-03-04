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

  plot, [0], [0], $
        xr=[0, 1.2], xstyle=1, xtitle='redshift', $
        yr=[0, 1], ystyle=1, ytitle=textoidl('vote fraction for question X, f_{x,z}'), $
        xtickv=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2], xticks=6, $
        ytickv=[0.0, 0.25, 0.50, 0.75, 1.0], yticks=5


; Some fake data

; Redshift array
  redshift_array_1=[0.030, 0.3+findgen(8)*0.1]

; Make fake vote fractions for 3 vote fraction levels, 3 SB levels, 1+8 redshift intervals, 7
; evolutionary corrections
  vote=dblarr(3, 3, 7, 9)

  loadct,26,/silent

; Loop over vote fraction levels levels
  for i=0L,3-1 do begin
     
                                ; Loop over SB levels
     for j=0L,3-1 do begin
        
                                ; Set original z~0 vote fraction levels
        vote(i,j,*,0)=0.28d*(i+1)+randomn(seed, 1)/30.
        offset=vote(i,j,*,0)
        offset=offset(0)


                                ; Loop over evolutionary corrections,
                                ; add random vote fractions
        for k=0L,7-1 do begin
           vote(i,j,k,1:8)= randomn(seed, 8)/([20d, replicate(25d, 8)])  + offset
           
                                ; for fake vote fractions: no
                                ; evolutionary correction at low-z
                                ; bin, so set them equal
           
           
                                ; Plot
           color=25+i*75 + k*3
           oplot, [redshift_array_1], [vote(i,j,k,*)], psym=sym(1), color=color
           oplot, [redshift_array_1], [vote(i,j,k,*)], linestyle=0, color=color
        endfor
     endfor
  endfor




; galaxies with an in-between vote
;oplot, [0.012], [0.50], psym=sym(1)
;oplot, [0.035], [0.56], psym=sym(1)
;oplot, [0.008], [0.47], psym=sym(1)


; galaxies with a high vote
;oplot, [0.016], [0.88], psym=sym(1)
;oplot, [0.030], [0.92], psym=sym(1)
;oplot, [0.022], [0.94], psym=sym(1)

; label
;xyouts, [0.23], [0.2], 'SDSS galaxies', /normal, orientation=90, charsize=1.2


;--------------------------------------------------------------------
; ## Clean up
  device, /close
  set_plot, 'x'
  cleanplot, /silent
;--------------------------------------------------------------------

end
