function colorbar_pos,pos,dx=dx,width=witdh
   if n_elements(dx)    ne 1 then dx    = (pos[2]-pos[0])*0.02d0
   if n_elements(width) ne 1 then width = (pos[2]-pos[0])*0.05d0
   pos   = [pos[2]+dx,pos[1],pos[2]+dx+width,pos[3]]
   if pos[2] gt 1d0 then pos[2] = 1d0
   return,pos
end
