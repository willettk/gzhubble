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




;--------------------------------------------------------------------
; ## Clean up
device, /close
set_plot, 'x'
cleanplot, /silent
;--------------------------------------------------------------------

end
