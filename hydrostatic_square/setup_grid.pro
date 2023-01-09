pro setup_grid,vx=vx,vy=vy,res=res,regular=regular,oname=oname

  if (n_elements(vx) eq 0) then vx = 142.3
  if (n_elements(vy) eq 0) then vy = -31.4

  if (n_elements(oname) eq 0) then oname="square"
  
  scaling = 1d0
  rot = 0 ; 0 if not, 1 if rotated
  
  box_length = 1d0

  P = 2.5
  P_box = P
  P_amb = P

  internal_box_length = 0.5d0

  gamma = 5.0/3.0 ; 1.4

  if (n_elements(res) eq 0) then res = 2*2 ; only even numbers!

  rho_box = 4 ; internal box (square) ; changed to make filling with particles easier
  rho_amb = 1 ; ambient (everything except square)

  fout=oname+".ic"

  if (not keyword_set(regular)) then begin
  Ng_box = 2L^res * 4L
  
  ; total number of particles
  Np_box =  2L^(2*res) * 16L
  Np_amb =  2L^(2*res) * 12L ; = sqrt(3 (area) / 3 (density)) * Ng_amb
  
  Np = Np_box + Np_amb
  print,Np

  pos = fltarr(3,Np)
  vel = fltarr(3,Np)
  mass = fltarr(Np)
  id = lindgen(Np)+1
  u = fltarr(Np)

  mass_box = rho_box*(scaling*internal_box_length^2)/Np_box
  mass_amb = mass_box ; = rho_amb*(box_length^2 - internal_box_length^2)/Np_amb 

  ; box particles
  for i=0L, Ng_box-1 do begin
     for j=0L, Ng_box-1 do begin

        ip = i*Ng_box + j
 
        pos(0,ip) = 0.25 + (i+0.5d0)*internal_box_length/Ng_box
        pos(1,ip) = 0.25 + (j+0.5d0)*internal_box_length/Ng_box
        vel(0,ip) = vx
        vel(1,ip) = vy
        u(ip)= P_box/(gamma-1)/rho_box
        mass(ip) = mass_box

     endfor
  endfor

  ; ambient particles
  x_positions = [0,0   ,0  ,0   ,0.25,0.25,0.5,0.5 ,0.75,0.75,0.75,0.75]
  y_positions = [0,0.25,0.5,0.75,0   ,0.75,0  ,0.75,0   ,0.25,0.5 ,0.75]
  Ng_amb = 2^res
  for k=0,11 do begin
     x_pos=x_positions[k]
     y_pos=y_positions[k]
     for i=0L, Ng_amb-1 do begin
        for j=0L, Ng_amb-1 do begin
  
           ip= 12*i*Ng_amb + 12*j + k + Np_box
        
           pos(0,ip)= 0.25*(i+0.5d0)*box_length/Ng_amb + x_pos 
           pos(1,ip)= 0.25*(j+0.5d0)*box_length/Ng_amb  + y_pos
           vel(0,ip) = vx
           vel(1,ip) = vy
           u(ip)= P_amb/(gamma-1)/rho_amb
           mass(ip) = mass_amb

        endfor
     endfor
  endfor 

  for i=0,Np-1 do begin
     pos(0,i) = pos(0,i)*scaling
     pos(1,i) = pos(1,i)*scaling
     pos(2,i) = pos(2,i)*scaling

     if (rot ne 0) then begin
        tmp1 = pos(0,i)
        pos(0,i) = pos(1,i)
        pos(1,i) = tmp1
     endif 
  endfor

endif else begin

   Np = 2L^(2*res)
   print,Np

   pos = fltarr(3,Np)
   vel = fltarr(3,Np)
   mass = fltarr(Np)
   id = lindgen(Np)+1
   u = fltarr(Np)
   
   mass_box = rho_box*scaling/Np
   mass_amb = rho_amb*scaling/Np 
   for i=0L, 2L^res - 1  do begin
      for j=0L, 2L^res - 1 do begin

         ip = i*2L^(res) + j
 
         pos(0,ip) = (i+0.5d0)/2L^res
         pos(1,ip) = (j+0.5d0)/2L^res
         vel(0,ip) = vx
         vel(1,ip) = vy
         if (pos(0,ip) gt 0.25) and (pos(0,ip) lt 0.75) and (pos(1,ip) gt 0.25) and (pos(1,ip) lt 0.75) then begin 
            u(ip)= P_box/(gamma-1)/rho_box
            mass(ip) = mass_box
         endif else begin
            u(ip)= P_amb/(gamma-1)/rho_amb
            mass(ip) = mass_amb
         endelse 
         
      endfor
   endfor

   for i=0,Np-1 do begin
      pos(0,i) = pos(0,i)*scaling
      pos(1,i) = pos(1,i)*scaling
      pos(2,i) = pos(2,i)*scaling
      
      if (rot ne 0) then begin
         tmp1 = pos(0,i)
         pos(0,i) = pos(1,i)
         pos(1,i) = tmp1
      endif 
   endfor
   
  endelse
     
  print, "mass=",minmax(mass)
  print, "u   =",minmax(u)
  print, "vel =",minmax(vel)
  
  npart=lonarr(6)
  massarr=dblarr(6)
  time=0.0D
  redshift=0.0D
  flag_sfr=0L
  flag_feedback=0L
  partTotal=lonarr(6)
  flag_cooling=0L
  num_files=1L
  BoxSize=1.0D
  Omega0=0.0D
  OmegaLambda=0.0D
  HubbleParam=0.7D
  flag_stellarage=0L            ; /*!< flags whether the file contains formation times of star particles */
  flag_metals=0L                ; /*!< flags whether the file contains metallicity values for gas and star particles */ 
  npartTotalHighWord=lonarr(6)  ;
  Labels=bytarr(2,15)
  bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4 - 6*4 - 2*4 - 4*8 - 2*4 - 6*4 - 2*15
  la=bytarr(bytesleft)
  
  h = { head , npart:npart,$
        massarr:massarr,$
        time:time,$
        redshift:redshift,$
        flag_sfr:flag_sfr,$
        flag_feedback:flag_feedback,$
        partTotal:partTotal,$
        flag_cooling:flag_cooling,$
        num_files:num_files,$
        BoxSize:BoxSize,$
        Omega0:Omega0,$
        OmegaLambda:OmegaLambda,$
        HubbleParam:HubbleParam,$
        flag_stellarage:flag_stellarage,$
        flag_metals:flag_metals,$
        npartTotalHighWord:npartTotalHighWord,$
        Labels:Labels,$
        la:la}
  
  h.npart(0)=Np
  h.partTotal(0)=Np

  write_head,fout,h
  add_block,fout,pos,'POS '
  add_block,fout,vel,'VEL '
  add_block,fout,id,'ID  '
  add_block,fout,mass,'MASS'
  add_block,fout,u,'U   '
  
end



