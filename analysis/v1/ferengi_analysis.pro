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
low1_redshift=[0.030, 0.3+findgen(8)*0.1]
low2_redshift=[0.030, 0.3+findgen(8)*0.1]
low3_redshift=[0.030, 0.3+findgen(8)*0.1]

low1_vote=[0.3+randomn(seed, 9)/([50d, replicate(30d, 8)])]
low2_vote=[0.3+randomn(seed, 9)/([50d, replicate(30d, 8)])]
low3_vote=[0.3+randomn(seed, 9)/([50d, replicate(30d, 8)])]

oplot, [low1_redshift], [low1_vote], psym=sym(1)
oplot, [low1_redshift], [low1_vote], linestyle=0

oplot, [low2_redshift], [low2_vote], psym=sym(1)
oplot, [low2_redshift], [low2_vote], linestyle=0

oplot, [low3_redshift], [low3_vote], psym=sym(1)
oplot, [low3_redshift], [low3_vote], linestyle=0
stop

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
