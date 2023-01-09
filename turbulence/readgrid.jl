# reads Sph2Grid output in HDF5 format

using HDF5

"""
returns Sph2Grid data
fname: input filename
block: blockname
fft: return complex FFT grid instead
     block "KVEC" gives k-vector grid
Pk: return binned bowerspectrum
    block "KPK" gives k-vector grid
head: Header strucutre, optional
scal: return length of a vector grid
"""
function readgrid(fname::String,block::String;fft::Bool=false,Pk::Bool=false,scalar::Bool=false)
  if cmp(block,"KSCA") == 0
    block="KVEC"
    scalar=true
  end

  file_id = h5open(fname,"r")

  head = read(file_id,"HEAD")

  if block == "HEAD"
    return head
  end

  if fft
    rdata=read(file_id, block*"/FFTGrid_real")
    idata=read(file_id, block*"/FFTGrid_imag")
    grid = rdata + im*idata

  elseif Pk
    grid = read(file_id, block*"/Pk")
  else
    grid = read(file_id, block*"/Grid")
  end

  return grid
  close(file_id)
  
end

