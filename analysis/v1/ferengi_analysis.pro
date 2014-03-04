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

; galaxies with a low vote
oplot, [0.010], [0.30], psym=sym(1)
oplot, [0.020], [0.27], psym=sym(1)
oplot, [0.015], [0.22], psym=sym(1)


; galaxies with an in-between vote
oplot, [0.012], [0.50], psym=sym(1)
oplot, [0.035], [0.56], psym=sym(1)
oplot, [0.008], [0.47], psym=sym(1)


; galaxies with a high vote
oplot, [0.016], [0.88], psym=sym(1)
oplot, [0.030], [0.92], psym=sym(1)
oplot, [0.022], [0.94], psym=sym(1)

; label
;xyouts, [0.23], [0.2], 'SDSS galaxies', /normal, orientation=90, charsize=1.2


;--------------------------------------------------------------------
; ## Clean up
device, /close
set_plot, 'x'
cleanplot, /silent
;--------------------------------------------------------------------

end
