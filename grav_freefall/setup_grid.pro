pro setup_grid

  scaling = 1d0
  
  box_length = 2d0
  radius = box_length * 0.5


; 20 per direction
  ; then cut out a sphere
  
  internal_box_length = 0.5d0

  gamma = 5.0/3.0

  P = 1e-6
  mass_sphere = 1.0
  
  fout="freefall.ic"

  Ng_box = 20
  Np = Ng_box^3
  print,Np

  Np_sphere = 0

  ; find out, how many of these lie inside the sphere
  for i=0, Ng_box-1 do begin
     for j=0, Ng_box-1 do begin
        for k=0, Ng_box-1 do begin
           x = (i+0.5)/Ng_box * box_length
           y = (j+0.5)/Ng_box * box_length
           z = (k+0.5)/Ng_box * box_length
           r2 = (x-radius)^2 + (y-radius)^2 + (z-radius)^2 ; radius = offset
           if r2 le radius^2 then Np_sphere = Np_sphere+1
        endfor
     endfor
  endfor

  print,Np_sphere
  rho = mass_sphere / 8.0 * Np / Np_sphere
  
  pos = fltarr(3,Np_sphere)
  vel = fltarr(3,Np_sphere)
  mass = fltarr(Np_sphere)
  id = lindgen(Np_sphere)+1
  u = fltarr(Np_sphere)
  rho_arr = fltarr(Np_sphere)

  ip = 0
  
  ; arrange particles
  for i=0L, Ng_box-1 do begin
     for j=0L, Ng_box-1 do begin
        for k=0L, Ng_box-1 do begin
           x = (i+0.5)/Ng_box * box_length
           y = (j+0.5)/Ng_box * box_length
           z = (k+0.5)/Ng_box * box_length
           r2 = (x-radius)^2 + (y-radius)^2 + (z-radius)^2 ; radius = offset
           if r2 le radius^2 then begin
              pos(0,ip) = x
              pos(1,ip) = y
              pos(2,ip) = z
              vel(0,ip) = 0
              vel(1,ip) = 0
              vel(2,ip) = 0
              u(ip)= P/(gamma-1)/rho
              mass(ip) = mass_sphere / Np_sphere
              rho_arr(ip) = rho

              ip = ip + 1
           endif
        endfor
     endfor
  endfor
     
  print, "mass=",minmax(mass)
  print, "u   =",minmax(u)
  print, "vel =",minmax(vel)
  
  npart=lonarr(6)
  massarr=dblarr(6)
  massarr(0) = mass_sphere / Np_sphere ; set the correct value also here
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
  flag_stellarage=0L            ;          /*!< flags whether the file contains formation times of star particles */
  flag_metals=0L                ;              /*!< flags whether the file contains metallicity values for gas and star particles */ 
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
  
  h.npart(0)=Np_sphere
  h.partTotal(0)=Np_sphere

  write_head,fout,h
  add_block,fout,pos,'POS '
  add_block,fout,vel,'VEL '
  add_block,fout,id,'ID  '
  add_block,fout,mass,'MASS'
  add_block,fout,u,'U   '
  add_block,fout,rho_arr,'RHO '
  
end



