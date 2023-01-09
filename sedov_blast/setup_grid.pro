pro setup_grid,single_particle=single_particle

  box_length = 1d0

  P = 1d-6
  P_amb = P

  gamma = 5.0/3.0

  if (not keyword_set(single_particle)) then begin
     Np = 64L                   ; particles per direction
     high_part0 = Np/2.-1
     high_part1 = Np/2
     n_dist = 8.0
     print,"distributed"
  endif else begin
     Np = 65L
     high_part0 = (Np-1)/2
     high_part1 = (Np-1)/2
     n_dist = 1.0
     print,"single particle"
  endelse
  
  Ng =  Np*Np*Np ; total number of particles

  rho = 1.0
  E = 10.0
  
  fout="sedov.ic"

  pos = fltarr(3,Ng)
  vel = fltarr(3,Ng)
  mass = fltarr(Ng)
  id = lindgen(Ng)+1
  u = fltarr(Ng)

  mass_particle = rho * box_length*box_length*box_length / Ng 

  ; setup grid
  for i=0L, Np-1 do begin
     for j=0L, Np-1 do begin
        for k=0L, Np-1 do begin
           ip = i*Np*Np + j*Np + k
 
           pos(0,ip) = (i+0.5d0)*box_length/Np
           pos(1,ip) = (j+0.5d0)*box_length/Np
           pos(2,ip) = (k+0.5d0)*box_length/Np
           vel(0,ip) = 0
           vel(1,ip) = 0
           vel(2,ip) = 0
           u(ip) = P_amb/(gamma-1)/rho
           mass(ip) = mass_particle
        endfor 
     endfor
  endfor

  ; initialize eigth particles with high internal energy (kernel)
  for i=high_part0, high_part1 do begin
     for j=high_part0, high_part1 do begin
        for k=high_part0, high_part1 do begin
           ip = i*Np*Np + j*Np + k
           print,pos(*,ip)
           
           u(ip) = E / n_dist / mass_particle
        endfor
     endfor
  endfor
  
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
  
  h.npart(0)=Ng
  h.partTotal(0)=Ng

  write_head,fout,h
  add_block,fout,pos,'POS '
  add_block,fout,vel,'VEL '
  add_block,fout,id,'ID  '
  add_block,fout,mass,'MASS'
  add_block,fout,u,'U   '
  
end



