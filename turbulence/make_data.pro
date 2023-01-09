; This file contains routines to produce a gadget snapshot 
; of a periodic box with a Kolmogorov velocity power spectrum
; from a glass file. The turbulent velocity contains the fraction
; X_turb of the thermal energy in the box. Temperature and density
; are set to values reasonable for a galaxy cluster atmosphere. 
; Note: This is not a real GADGET snapshot as hsml is set to a 
; reasonable constant. Upon running GADGET is going to recalculate it.
;
; Some routines from Klaus' library are needed:
; * write_head.pro
; * add_block.pro
;
; Note: you have to compile this twice
; originally produced by P. Foerster

; make initial conditions for a turbulent box
pro make_box, npart, debug=debug,X_turb=X_turb

	if not keyword_set(npart) then $
		npart = 128L^3

        if not keyword_set(X_turb) then $
		X_turb = 0.3

	npart = ulong(npart)

	print, "Npart = ", npart

; Gadget units & chemistry
    m_unit = 1.989d43       ; [10^10 Msol]
    l_unit = 3.085678d21    ; [kpc]
    v_unit = 100000D        ; [km/s]
    
    H_frac = 0.76
	umol = 4.0/(5.0*H_frac+3.0)
    mp = 1.6726231e-24      ; proton mass cgs
    k_boltz = 1.3806580e-16 ; [cgs]

; input values
    fout = './snap_000'     ; output filename
    ;X_turb = 0.3            ; E_turb / E_therm
    
    boxsize = 3000D         ; [kpc]
    
    T = 1d7                 ; [K]
    rho = 1d-27/(m_unit/l_unit^3.)  ; [GADGET]
    mass = rho * boxsize^3
    
; make positions
    pos = make_positions_hcp(boxsize, boxsize, boxsize, npart)

; make data structures
    head = make_head()
	vel	 = make_array(3,npart,/float)
	id	 = make_array(npart,/uint)
	u	 = make_array(npart,/float)
	hsml = make_array(npart,/float)
	dens = make_array(npart,/float)

    head.npart = [1,0,0,0,0,0] * npart
	head.massarr = [1,0,0,0,0,0] * mass/npart
	head.time = 1
	head.redshift = 1./head.time -1
	head.flag_sfr = 0
	head.flag_feedback = 0
	head.parttotal = head.npart
	head.flag_cooling = 0
	head.num_files = 1
	head.boxsize = boxsize
	head.omega0 = 0.3
	head.omegalambda = 0.7
	head.hubbleparam = 0.7

; IDs
    id = ulindgen(npart)+1

; internal energy
    u[*] = u2t(T,/inv)

; thermal energy
    mass_cgs = mass * m_unit
	Etherm = mass_cgs/umol/mp * 3./2. *k_boltz * T

    print, 'Thermal Energy in Box [cgs] = ',Etherm

; hsml
    hsml[*] = 75.51 ; what GADGET finds on average for this glass

; make velocities (the hard part) 
	Etherm /= m_unit *v_unit^2         ; cgs to gadget
    v2 = npart*2.0*Etherm/mass*X_turb      ; total v^2 of particles
    print, "Total amplitude of (squared) velocities = ",v2 
    
    ngrid = long(npart^(1./3.))

    vgrid = make_vel(ngrid, boxsize, v2, debug=debug)
    
    ; Ncells = Npart 
    cellsize = boxsize / ngrid

    vel[0,*] = idlNGP(pos/cellsize, vgrid[0,*,*,*])
    vel[1,*] = idlNGP(pos/cellsize, vgrid[1,*,*,*])
    vel[2,*] = idlNGP(pos/cellsize, vgrid[2,*,*,*])

; density
    dens[*] = rho

; output
    print, 'Writing : '+fout

    write_head, fout, head
	add_block, fout, float(pos), 'POS'
	add_block, fout, float(vel), 'VEL'
	add_block, fout, ulong(id), 'ID'
	add_block, fout, float(u), 'U'
	add_block, fout, float(hsml), 'HSML'
	add_block, fout, float(dens), 'RHO'

end

; make velocity grid (this is where the magic happens :-)
function make_vel, ngrid, boxsize, amp, debug=debug
 
    seed = 14041981

	kmin = 2*!pi/(boxsize)	    ; box mode
	kmax = !pi*ngrid/boxsize    ; Nyquist mode

	vel = make_array(3,ngrid, ngrid, ngrid,/float, val=0)	
    kmag = make_array(ngrid,ngrid,ngrid,/float, val=0)
    cdata = make_array(3,ngrid,ngrid,ngrid,/complex, val=0)
    cdata_rl = make_array(ngrid,ngrid,ngrid, val=0,/double)
    cdata_im = make_array(ngrid,ngrid,ngrid, val=0,/double)
    iconj = 0
    jconj = 0
    kconj = 0

	for axes=0,2 do begin
		for i=0,ngrid-1 do $
		for j=0,ngrid-1 do $
		for k=0,ngrid/2 do begin
; Generate k value	first (thank you Volker)
            ; Define conjugated indizes of the grid
			if i ne 0 then iconj = ngrid - i $
					  else iconj = 0
			if j ne 0 then jconj = ngrid - j $
					  else jconj = 0
			if k ne 0 then kconj = ngrid - k $
					  else kconj = 0

            ; Define grid
			if i LE ngrid/2. then kx = i * kmin $
							 else kx = -iconj * kmin
	
			if j LE ngrid/2. then ky = j * kmin $
							 else ky = -jconj * kmin
	
			if k LE ngrid/2. then kz = k * kmin $
						     else kz = -kconj * kmin

			kmag[i,j,k]  = sqrt(kx^2 + ky^2 + kz^2)
	
 
			if kmag[i,j,k] GT kmax then $
				continue    ; Only do a sphere in k space
	
            if i+j+k eq 0 then $
                continue    ; no DC current
           ; if (i eq ngrid/2) or (j eq ngrid/2) or (k eq ngrid/2) then $
           ;     continue    ; no DC current

            if i gt ngrid/2 then $ 
                continue    ; these are done via symmetry

; Power spectrum P(k)
		    Pk = kmin* kolmog_3D(kmag[i,j,k],kmax,kmin) 

; Generate normal distributed random numbers with dispersion Pk
; using Box Mueller method
		 	A =  sqrt( -alog(randomu(seed,/double)) * Pk ) 
			phase = 2.*!pi*randomu(seed,/double)

; Cutting off all except the ~70 largest modes (Bauer&Springel2012)
                        IF kmag[i,j,k] LT (6.25/4)*kmin THEN $
                           A = 0
                        IF kmag[i,j,k] GT (12.57/4)*kmin THEN $
                           A = 0
	
; Set power so we get a real vel after inverse FFT
            if i gt 0 then begin    ; grid is hermitian in i>ngrid/2
                cdata_rl[i,j,k] = A * cos(phase)
		    	cdata_im[i,j,k] = A * sin(phase)

                cdata_rl[iconj,jconj,kconj] = cdata_rl[i,j,k]
    			cdata_im[iconj,jconj,kconj] = -1*cdata_im[i,j,k]
            end else begin  ; i = 0 needs special treatment
                if j eq 0 then begin ; first row
                    
                    if k gt ngrid/2. then $
                        continue
                
                    cdata_rl[i,j,k] = A * cos(phase)
		    	    cdata_im[i,j,k] = A * sin(phase)
                
                    cdata_rl[i,j,kconj] = cdata_rl[i,j,k]
    			    cdata_im[i,j,kconj] = -1*cdata_im[i,j,k]
                end else begin  ; j != 0 here
                    if j gt ngrid/2. then $     ; rest of the plane
                        continue
                
                    cdata_rl[i,j,k] = A * cos(phase)
		    	    cdata_im[i,j,k] = A * sin(phase)

                    cdata_rl[i,jconj,kconj] = cdata_rl[i,j,k]
    		    	cdata_im[i,jconj,kconj] = -1*cdata_im[i,j,k]
                end
            end
		end

        cdata[axes,*,*,*] = COMPLEX(cdata_rl,cdata_im)
        data = reform(FFT(cdata[0,*,*,*], /inverse, /double))/Ngrid^3

; check if we got the symmetries correct
        for i=0,ngrid^3-1 do $
            if abs(imaginary(data[i])) gt 1e-7*abs(real_part(data[i])) then $ 
                stop

; set velocities
		vel[axes,*,*,*] =  real_part(data)

	end

; norm to total amplitude because kolmog_3D is wrong
        norm = sqrt(total(vel[0,*,*,*]^2 +vel[1,*,*,*]^2 +vel[2,*,*,*]^2 ))
        
        vel *=  sqrt(amp) / norm
        cdata *=  sqrt(amp) / norm  ; Parseval's theorem -> ngrid^3

; testing
    if keyword_set(debug) then begin
        ; total energy of data 
        v2_k = total(abs(FFT(vel[0,*,*,*]))^2. $
                    +abs(FFT(vel[1,*,*,*]))^2. $
                    +abs(FFT(vel[2,*,*,*]))^2.) * ngrid^3
	    print, 'v^2  in k space : '+strn(v2_k)

        ; total energy in real space
	    v2 = total(vel[0,*,*,*]^2 +vel[1,*,*,*]^2 +vel[2,*,*,*]^2 )
    	print, 'v^2 in real space: '+strn(v2)
	    print, 'Requested v^2  : '+strn(amp)

        ; 3D Spectrum from data
        lcData = cData[0,*,*,*]^2+cData[1,*,*,*]^2+cData[2,*,*,*]^2

        plot, kmag,sqrt(lcdata *4*!pi * kmin^2 )  $
            ,xrange=[2e-3, 0.2],psym=3,yrange=[1e-8,1e3],/ylog,/xlog $ 
            ,ytitle='sqrt(P(k) 4!7p!X k!Dmin!N!U2!N) [km/s]',xstyle=1 $
            ,xtitle='k=2!7p!X/L [1/kpc]'
    	oplot, [1.,1]*kmin, [1e-10,1e10],col=16711680
	    oplot, [1.,1]*kmax, [1e-10,1e10],col=16711680

        ; 3D spectrum analytically
	    k = kmin + findgen(10000)/9999. * (kmax-kmin)
    	Pk =  kolmog_3D(k, kmax, kmin) * 8
    	oplot, k, sqrt(Pk * 4 * !pi * kmin^2), col=65280

        ; Spectrum from vel, this leaks into large k, because of numerical error
        kData = make_array(3,ngrid,ngrid,ngrid)
    	kData[0,*,*,*] = abs(FFT(vel[0,*,*,*],/double)) * ngrid^3
    	kData[1,*,*,*] = abs(FFT(vel[1,*,*,*],/double)) * ngrid^3
	    kData[2,*,*,*] = abs(FFT(vel[2,*,*,*],/double)) * ngrid^3
        
        lkData = kData[0,*,*,*]^2 + kData[1,*,*,*]^2 + kData[2,*,*,*]^2
    
        oplot, kmag, sqrt(lkdata * 4 * !pi * kmin^2)  ,col=16711680,psym=3

        stop
    end 

	return, vel 
end

function kolmog_3D, k, kmax, kmin 
; here we require 1 = int^kmax_kmin dk P_0 * 4 pi k^2 k^(-11/3)
    norm = (6*!pi * ( kmin^(-2./3.) - kmax^(-2./3.) ) )^(-1) 
	return, norm * k^(-11.0/3.0) 
end

; sample Grid at Pos via NGP the IDL way
function idlNGP, pos, ingrid

    ngrid = n_elements(ingrid)^(1./3.)
    
    grid = reform(ingrid,ngrid,ngrid,ngrid )

    npos = n_elements(pos)/3

    u = reform(pos[0,*],npos)
    v = reform(pos[1,*],npos)
    w = reform(pos[2,*],npos)

    i = floor(u)
    j = floor(v)
    k = floor(w)

    print,n_elements(i)
    print,n_elements(j)
    return, grid[i,j,k]
end

; make a gadget header
function make_head
	head = {        npart		        :	lonarr(6),	$
			massarr			:	dblarr(6),	$
			time			: 	double(1),	$
			redshift		:	double(1),	$
			flag_sfr		:	long(0),	$
                        flag_feedback	        :	long(0), 	$
			parttotal		: 	lonarr(6),	$
			flag_cooling	        :	long(0),	$
			num_files		:	long(1),	$
			boxsize			:	double(1),	$
			omega0			:	double(0.3),$
			omegalambda		:	double(0.7),$
			hubbleparam		:	double(0.7),$
                        flag_stellarage         :       long(0),$
                        flag_metals             :       long(0), $
                        npartTotalHighWord      :       lonarr(6), $
		        Labels                  :       bytarr(2,15),$
                        la                      :       bytarr(256-6*4 - 6*8 - 8 - 8 - 2*4 - 6*4 - 2*4 - 4*8 - 2*4 - 6*4 - 2*15)}

	return, head
	
end

; convert internal energy to temperature for GADGET
FUNCTION u2t, u,  rad=rad, inv=inv, xH=xH, uvel=uvel, gamma=gamma
	
    if not keyword_set(gamma) then $
        gamma = 5./3.

	IF NOT keyword_set(xH) THEN $
		xH = 0.76

	IF NOT keyword_set(uvel) THEN $
		uvel=1e5

	bk=1.380658d-16   	    ;k_boltzmann in cgs
	prtn=1.672623d-24		;m_proton in g 

	yhelium = ( 1. - xH ) / ( 4 * xH ) 

	mean_mol_weight = (1. + 4. * yhelium) / (1. + 3. * yhelium + 1) 

	IF keyword_set(inv) THEN BEGIN
		T = U
		u = T /( (gamma-1) * uvel^2 * prtn * mean_mol_weight / bk) 
		return, u
	END ELSE BEGIN
		T = u * (gamma-1) * uvel^2 * prtn * mean_mol_weight / bk
		return, T
	END
END

; optimal hcp particle distribution in a periodic 
; box, with ntot being a rough upper bound of particles to
; distribute. return pos and ntot
; ntot values like X^3 recommended for cubic boxes
function make_positions_hcp, Lx, Ly, Lz, ntot, debug=debug

    if n_params() lt 3 then begin
        print, 'Usage : pos = make_positions_hcp(Lx, Ly, Lz, ntot, debug=debug)'
        print, '  Lx    : Boxsize in x direction'
        print, '  Ly    : Boxsize in y direction'
        print, '  Lz    : Boxsize in z direction'
        print, '  ntot  : Max total number of particles'
        print, '  debug : Show periodicity diagnostics'
        return, -1
    end

    ; ntot is better divisible by two
    if ntot mod 2 ne 0 then $
        ntot--

    ; find by combining spacings with Ntot
    r = (sqrt(2.0D)*Lx*Ly*Lz/8.0D /ntot)^(1.0D/3D)

    ; spacings
    dx = 2.0D*r
    dy = sqrt(3.0D)*r
    dz = sqrt(6.0D)*2.0D/3.0D*r
   
    ; particle numbers
	np = make_array(3, val=0D)
	np[*] = long(ntot^(1./3.))
	
    ; enforce periodicity
    dx += (lx-dx*np[0])/np[0]
    dy += (ly-dy*np[1])/np[1]
    dz += (lz-dz*np[2])/np[2]
       
    ; diagnostics
    print, 'Particle numbers : '
    print, '     Nx = '+strn(np[0])
    print, '     Ny = '+strn(np[1])
    print, '     Nz = '+strn(np[2])
    print, ' Total  = '+strn(np[0]*np[1]*np[2])
    print, ' Wanted = '+strn(ntot)
    print, ' Delta  = '+strn(double(ntot)-np[0]*np[1]*np[2])

    ntot = np[0]*np[1]*np[2]

    ; particle positions
    x = make_array(np[0],np[1],np[2],/double)
    y = make_array(np[0],np[1],np[2],/double)
    z = make_array(np[0],np[1],np[2],/double)

    idxarr = double(lindgen(np[0]))

    ; A(0) plane
    ; 0st row
    x[*,0,0] = r + idxarr*dx
    y[*,0,0] = r
    z[*,0,0] = r
    
    ; 1st row
    x[*,1,0] = idxarr*dx
    y[*,1,0] = r + dy
    z[*,1,0] = r

    ; A-plane
    ; even rows
    for i=2L,np[1]-1,2 do begin
        x[*,i,0] = x[*,0,0]
        y[*,i,0] = y[*,0,0] + i*dy
        z[*,i,0] = z[*,0,0]
    end

    ; odd rows
    for i=3L,np[1]-1,2 do begin
        x[*,i,0] = x[*,1,0]
        y[*,i,0] = y[*,1,0] + (i-1)*dy
        z[*,i,0] = z[*,1,0]
    end

    ;B(1)-plane
    ; 0rst row
    x[*,0,1] = r + idxarr*dx 
    y[*,0,1] = 0 
    z[*,0,1] = r + dz

    ; 1st row
    x[*,1,1] = idxarr*dx 
    y[*,1,1] = dy
    z[*,1,1] = r + dz

    ; B-plane
    ; even rows
    for i=2L, np[1]-1, 2 do begin
        x[*,i,1] = x[*,0,1]
        y[*,i,1] = y[*,0,1] + i*dy
        z[*,i,1] = z[*,0,1]
    end
   
    ; odd rows
    for i=3L, np[1]-1, 2 do begin
        x[*,i,1] = x[*,1,1]
        y[*,i,1] = y[*,1,1] + (i-1)*dy
        z[*,i,1] = z[*,1,1]
    end
    
    ; all planes
    ; even planes
    for i=2L, np[2]-1, 2 do begin
        x[*,*,i] = x[*,*,0]
        y[*,*,i] = y[*,*,0]
        z[*,*,i] = z[*,*,0] + i*dz
    end
    
    ; odd planes
    for i=3L, np[2]-1, 2 do begin
        x[*,*,i] = x[*,*,1]
        y[*,*,i] = y[*,*,1]
        z[*,*,i] = z[*,*,1] + (i-1)*dz
    end
    
    ; make particle arrays
    pos = make_array(3, ntot,/double)

    bin = ulong(0)
    for i=0,np[0]-1 do $
    for j=0,np[1]-1 do $
    for k=0,np[2]-1 do begin
        pos[0,bin] = x[i,j,k] + dx/4d0
        pos[1,bin] = y[i,j,k] + dy/4d0
        pos[2,bin] = z[i,j,k] + dz/4d0
        bin++
    end

    ; randomize
    seed = 14041981
    for k=0, 2*ntot do begin
        i = round(randomu(seed)*(ntot-1))
        j = round(randomu(seed)*(ntot-1))

        tmp = pos[*,i]
        pos[*,i] = pos[*,j]
        pos[*,j] = tmp
    end

    if keyword_set(debug) then begin
        old_pmulti = !p.multi
        !p.multi[1:2] = [3,2]

        ; check spacings
        plot, x, y, psym=4, /iso, xrange=[0, Lx], yrange=[0, Ly], $
            xtitle='x', ytitle='y'
        oplot, x[0,0,0]+[0,dx], y[0,0,0]*[1.,1.], col=color(1)

        plot, x, z, psym=4, /iso, xrange=[0, Lx], yrange=[0, Lx], $
            xtitle='x', ytitle='z'
        oplot, x[0,0,0]*[1,1], z[0,0,0]+[0,dz], col=color(1)

        plot, y, z, psym=4, /iso, xrange=[0, Ly], yrange=[0, Lz], $
            xtitle='y', ytitle='z'
        oplot, y[0,0,0]+[0,dy], z[0,0,0]*[1.,1.], col=color(1)

        ; check periodicity
        bad = where(x gt Lx or x lt 0, nbadx)
        bad = where(y gt Ly or y lt 0, nbady)
        bad = where(z gt Lz or z lt 0, nbadz)
        print, 'Points outside the box  :'+strn(nbadx+nbady+nbadz)

        err = make_array(3, /double)
        err[0] = Lx - dx*np[0]
        err[1] = Ly - dy*np[1]
        err[2] = Lz - dz*np[2]

        print, "Error in sampling periodicity : "
        print, "    "+strn(err[0],len=7)+" particle spacings in x"
        print, "    "+strn(err[1],len=7)+" particle spacings in y"
        print, "    "+strn(err[2],len=7)+" particle spacings in z"

        x[*,*,*] += 2*dx
        bad = where(x gt Lx)
        x[bad] -= Lx

        y[*,*,*] += 2*dy
        bad = where(y gt Ly)
        y[bad] -= Ly
        
        z[*,*,*] += 2*dz
        bad = where(z gt Lz)
        z[bad] -= Lz

        plot, x, y, psym=4, /iso, xrange=[0, Lx], yrange=[0, Ly], $
            xtitle='x', ytitle='y'
        oplot, 2*dx * [1,1] , [0.,Ly], col=color(1)
        oplot, [0,Lx] , 2*dy*[1,1], col=color(1)

        plot, x, z, psym=4, /iso, xrange=[0, Lx], yrange=[0, Lx], $
            xtitle='x', ytitle='z'
        oplot, [0,Lx] , 2*dz*[1,1], col=color(1)
        oplot, 2*dx * [1,1] , [0.,Lz], col=color(1)

        plot, y, z, psym=4, /iso, xrange=[0, Ly], yrange=[0, Lz], $
            xtitle='y', ytitle='z'
        oplot, 2*dy * [1,1] , [0.,Lz], col=color(1)
        oplot, [0,Ly] , 2*dz*[1,1], col=color(1)

        stop

        !p.multi = old_pmulti

        return, -1
    end

    return, float(pos)
end
