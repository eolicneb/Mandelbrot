

subroutine color(mu, pixel)
    real, intent(in) :: mu
    integer*2, dimension(1:3), intent(out) :: pixel
    real :: corte = 256.
    pixel(1) = int(corte/2 -1 + corte/2*sin(mu*2*3.1416/corte)) !mu*(mu-corte/2.)*(corte-mu)
    pixel(2) = int(corte/2 -1 + corte/2*sin((mu-corte/3.)*2*3.1416/corte)) !(mu+corte/3.)*(mu-corte/6.)*(corte*2./3.-mu)
    pixel(3) = int(corte/2 -1 + corte/2*sin((mu+corte/3.)*2*3.1416/corte)) !(mu-corte/3.)*(mu-5.*corte/6.)*(corte*4./3.-mu)
    !pixel(:) = int(mu)
end subroutine color

program mandelbrot
use bmp
implicit none
! gfortran -O2 mandel.f95 -o mandel
! ./mandel |./DynamicLattice -norange -nx 900 -ny 650 -pix 1
integer, parameter:: real12=selected_real_kind(16,16)
integer, parameter:: ancho = 900, alto = 650, margen = 20470, escape = 30
real(8), parameter:: Scent = -0.80055004, Tcent = 0.170759856, Diag = .1 !.000000004
integer:: i, j, n, escape2
integer*2,dimension(1:2048,1:2048,1:3) :: out=0
real :: vlr
real(8):: SMax, SMin, TMax, TMin, s, t, o, p, Zr, Zi, z
integer*2, dimension(1:3) :: pix
character*20 :: archivo = "mandel026.bmp"

SMax = Scent+ancho*Diag/(ancho**2+alto**2)**0.5
SMin = Scent-ancho*Diag/(ancho**2+alto**2)**0.5
TMax = Tcent+alto*Diag/(ancho**2+alto**2)**0.5
TMin = Tcent-alto*Diag/(ancho**2+alto**2)**0.5
o = (SMax-SMin)/(ancho-1)
p = (TMax-TMin)/(alto-1)
escape2=escape**2

do i=0,ancho-1
  do j=0,alto-1
    s = SMin+o*i; t = TMin+p*j
    Zr = 0.; Zi = 0.; z=0.
    n=0
    do while (Zr**2+Zi**2<escape2 .and. n<margen)
      n=n+1
      z = Zr
      Zr = (z**2-Zi**2)+s
      Zi = 2*z*Zi+t
      z=Zr**2+Zi**2
    enddo
    !z = exp((Zr**2 + Zi**2)**0.5)
    if (n<margen) then
      !mu = m + frac = n + 1 - log (log  |Z(n)|) / log 2
      vlr=n+1-log(log(sqrt(z)))/log(2.)
      !out(j+1,i+1,1)=vlr
      !out(j+1,i+1,3)=255-vlr
      call color(vlr,pix)
      out(j+1,i+1,:)=pix
    endif
  enddo
enddo
write(*,*) "Guardando en "//archivo
call saveBMP(archivo,out,ancho,alto,1,1,1)

end program mandelbrot
