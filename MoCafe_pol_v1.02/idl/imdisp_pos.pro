@imdisp.pro
function imdisp_pos,image,margin=margin,aspect=aspect,position=position
;+
; Usage:
;    pos = imdisp_pos(image,margin=,aspect=)
;
; Obtain position values used when image is plotted with imdisp.pro
;
; Written by K.-I. Seon/KASI
;
; Example 1):
;    !p.multi=[0,1,2]
;    image = findgen(100,100)
;    x = findgen(10)
;    y = x^2
;    plot,x,y,position=imdisp_pos(image)
;    imdisp,image,/eras,/axis
;
; Example 2):
;    !p.multi=[0,1,2]
;    image = findgen(100,100)
;    imdisp,image,/eras,/axis,position=pos
;    x = findgen(10)
;    y = x^2
;    plot,x,y,position=imdisp_pos(image,position=pos)
;-
  plot,[0],/nodata,/noerase,xstyle=4,ystyle=4,xmargin=[0,0],ymargin=[0,0]
  ;plot,[0],/nodata,/noerase,xstyle=4,ystyle=4
  if n_elements(position) eq 4 then begin
     pos = position
  endif else begin
     pos = [!x.window[0], !y.window[0], !x.window[1], !y.window[1]]
  endelse
  if n_elements(margin) eq 0 then margin = 0.025
  imdisp_imsize,image,x0,y0,xsize,ysize,position=pos,aspect=aspect,margin=margin
  return,pos
end
