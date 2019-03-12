!--------------------------------DEFINE FUNCTIONS------------------------------
function sgn(i,nK)

	integer :: i, sgn
	if (i <= nK/2 ) then     !!!!need change!!!#####################
		sgn = -1
	else
		sgn = 1
	end if
	return
end function sgn





!--------------------------------DEFINE FUNCTIONS------------------------------



program legleg

implicit none


integer, parameter :: extra = 0
!nK = (1 + 2 * extra)*4  !extra points
integer, parameter :: nK = (1 + 2 * extra)*4

real*8 kf1, kf2, vk1, dk, kspan

integer ::k1, k2, k3, k4, k5, k6
integer i1, i2, i3, i4, iq, ival, iqmax, i
integer, dimension(nK, nK, nK) ::leg4 = 0	
real*8, dimension(nK) :: Ek = 0, vK, Kvec   ! fermi energy, fermi velocity, fermi momentum
real*8, dimension(nK, nK) :: xpp, xph 		
real*8, dimension(nK, nK, nK) :: dGa = 0, Ga = 0
real*8 kv, Qvec, radius
real*8 b11r,b11s,b12r,b12s,f12r,f12s,u11r,u12r,u12s

real*8 U, Vnn, Jnn, tt
integer marki, markj, markk

real*8 :: twopi = asin(1.d0) * 4, pi = asin(1.d0) * 2
real*8 :: cutoff = 1.d-8, wUV = 100, divergence = 100
real*8 b, La, Lc, dL
real*8, external :: lindhard_1d
integer :: sgn

logical bandinfo

!analysis define

integer, parameter :: nq = 2
real*8, dimension(nq) :: qv
real*8 V(nK, nK), eval(nK)
integer, external :: findpatch
real*8 qsc, qsdw, qcdw, Vsc, Vsdw, Vcdw


qv = (/0.d0,pi/)

ival = 1;

!band info
bandinfo = .true.

if (bandinfo .eqv. .true. ) then

kf1 = -0.7 * pi; vk1 = 1.93649; kf2 = -0.3 * pi
kspan = .2 * pi

do i = 1, nK/4
	Kvec(i) = kf1 - kspan + (2*i-1) * kspan / (nK/4)
	Ek(i) = - ( Kvec(i)  - kf1 ) * vk1
	vK(i) = - vk1
	
end do

    !band 1, R movers
do i = nK/4 + 1, nK/2;  i1 = nK/2 - i + 1
	Kvec(i) = - Kvec(i1)
	Ek(i) = Ek(i1)
	vK(i) = - vK(i1)
end do

!band 2 is a copy of band 1, except for the shift of the real fermi points
dk = kf2 - kf1
Kvec(nK/2+1 : nK/2+nK/4) = Kvec(1:nK/4) + dk
Kvec(nK/2+nK/4+1 : nK) = Kvec(nK/4+1:nK/2) - dk
vK(nK/2+1 : nK) = vK(1:nK/2)
Ek(nK/2+1 : nK) = Ek(1:nK/2)

else

!---------------------------------------------------------
kf1 = -0.7 * pi; kf2 = -0.3 * pi
tt = 2. * cos(kf2)
! ek2 = -2 cos(k) + tt
! ek1 = -2 cos(k) - tt
kspan = .1 * pi

print*, tt, vk1

do i = 1, nK/4
	Kvec(i) = kf1 - kspan + (2*i-1) * kspan / (nK/4)
	Ek(i) = -2. * cos(Kvec(i)) - tt
	vK(i) = 2. * sin(Kvec(i))
	
end do

    !band 1, R movers
do i = nK/4 + 1, nK/2;  i1 = nK/2 - i + 1
	Kvec(i) = - Kvec(i1)
	Ek(i) = Ek(i1)
	vK(i) = - vK(i1)
end do

!band 2 is a copy of band 1, except for the shift of the real fermi points
dk = kf2 - kf1
Kvec(nK/2+1 : nK/2+nK/4) = Kvec(1:nK/4) + dk
Kvec(nK/2+nK/4+1 : nK) = Kvec(nK/4+1:nK/2) - dk
vK(nK/2+1 : nK) = vK(1:nK/2)
Ek(nK/2+1 : nK) = Ek(1:nK/2)
!-------------------------------------------------------------------


end if
!open(10, file='gaa.dat')
!do i = 1, nK
!	write(10,100) Kvec(i), Ek(i), vK(i)
!end do
!---------------------------------------------------------------------


!Kvec = (/-0.942478, -0.985335, -0.899621, -1.02819, -0.856764, -1.07105, -0.813906, 0.942478, 0.899621, 0.985335, 0.856764, 1.02819, 0.813906, 1.07105, -2.19911, -2.24197, -2.15626, -2.28483, -2.1134, -2.32769, -2.07054, 2.19911, 2.15626, 2.24197, 2.1134, 2.28483, 2.07054, 2.32769/)
!vK = (/-1.61803, -1.66691, -1.56618, -1.71273, -1.51145, -1.75541, -1.45395, 1.61803, 1.56618, 1.66691, 1.51145, 1.71273, 1.45395, 1.75541, -1.61803, -1.56618, -1.66691, -1.51145, -1.71273, -1.45395, -1.75541, 1.61803, 1.66691, 1.56618, 1.71273, 1.51145, 1.75541, 1.45395/)
!Ek = (/0., 0.0704025, -0.0682436, 0.142835, -0.134203, 0.217163, -0.197757, 0., -0.0682436, 0.0704025, -0.134203, 0.142835, -0.197757, 0.217163, 0., 0.0682436, -0.0704025, 0.134203, -0.142835, 0.197757, -0.217163, 0., -0.0704025, 0.0682436, -0.142835, 0.134203, -0.217163, 0.197757/)
!---------------------------------------------------------------------









!dk = (kf2 - kf1)/(2. * extra + 1.)
!dk = 0.3/(2. * extra + 1.)

radius = kspan/(nK/4)

!momentum conservation
Qvec = twopi
leg4 = 0
do k3 = 1, nK; do k2 = 1, nK; do k1 = 1, nK
	kv = mod (Kvec(k1) + Kvec(k2) - Kvec(k3), Qvec)
	if (kv > pi) kv = kv - twopi; if(kv < -pi) kv = kv + twopi
	!do k4 = 1, nK; if( abs(mod(kv - Kvec(k4),twopi)) >= dk / 2.) cycle
	do k4 = 1, nK; if( abs(mod(kv - Kvec(k4),twopi)) > 1.e-4) cycle
		leg4(k1,k2,k3) = k4; end do
end do; end do; end do


!assign vertex
U = 0.25/2;  Vnn = 0.10/2;  Jnn = 4 * ( U + Vnn )           !SC + SDW degenerate manifold, or SZH SO(5) manifold 
!U = 0.3;  Vnn = 0;  Jnn = 0.4                              !SC + CDW  degenerate manifold
!U = 0.3/2;  Vnn = -0.1/2;  Jnn = 0.8/2                     !SDW + CDW degenerate manifold
Vnn = Vnn - Jnn/4
Ga = 0



do i3 = 1, nK;  do i2 = 1, nK;  do i1 = 1, nK
	i4 = leg4(i1,i2,i3); if(i4 <= 0)cycle
	
	Ga(i1, i2, i3) = Ga(i1, i2, i3) - U * ( 1 + sgn(i1,nK) * sgn(i2,nK) * sgn(i3,nK) * sgn(i4,nK) ) / 4.
	Ga(i1, i2, i3) = Ga(i1, i2, i3) - Vnn * ( sgn(i2,nK) * sgn(i3,nK) + sgn(i1,nK) * sgn(i4,nK) ) / 4.        !coulomb interaction on rung
	Ga(i1, i2, i3) = Ga(i1, i2, i3) + 0.5 * Jnn * ( sgn(i1,nK) * sgn(i3,nK) + sgn(i2,nK) * sgn(i4,nK) ) / 4.  !antiferro spin exchange on rung (written as hop * hop which introduces a density-density interaction)

end do; end do;  end do


!assign vertex
!do k3 = 1, nK; do k2 = 1, nK; do k1 = 1, nK
!	k4 = leg4(k1, k2, k3); if(k4 == 0) cycle
!	Ga(k1, k2, k3) = -U
!end do; end do; end do








open(10, file = 'flowcurrent.dat')
open(20, file='flowVxxxx.dat')

La = wUV;  b = exp( log(wUV/cutoff) / 1200 )
do while (La > cutoff);  print*, La

	Lc = La * (1 + 1/b) / 2;  dL = La - La/b;  La = La / b

	!calculate xpp and xph
	do k5 = 1, nK;  do k6 = 1, nK
		xpp(k5, k6) = - lindhard_1d(Lc, Ek(k5), vK(k5), -Ek(k6), vK(k6), radius)
		xph(k5, k6) = lindhard_1d(Lc, Ek(k5), vK(k5), Ek(k6), vK(k6), radius)
    end do;  end do


!!!!!Feynamn diagram
!----------------------------------------------------
dGa = 0
do k3 = 1, nK; do k2 = 1, nK; do k1 = 1, nK
	k4 = leg4(k1,k2,k3); if(k4 == 0) cycle

	
	do k5 = 1, nK
		!PP channel
		k6 = leg4(k1, k2, k5)
		if( k6 /= 0) dGa(k1, k2, k3) = dGa(k1, k2, k3) + Ga(k1, k2, k5) * xpp(k5, k6) * Ga(k6, k5, k3)

		!PH,c channel
		k6 = leg4(k1, k5, k3)
        if(k6 /= 0) dGa(k1, k2, k3) = dGa(k1, k2, k3) + Ga(k1, k5, k3) * xph(k5, k6) * Ga(k6, k2, k5)

        !PH,d channel
        k6 = leg4(k5, k2, k3)
        if(k6 /= 0) dGa(k1, k2, k3) = dGa(k1, k2, k3) - 2 * Ga(k1, k6, k5) * xph(k5, k6) * Ga(k5, k2, k3)

        k6 = leg4(k5, k2, k3)
        if(k6 /= 0) dGa(k1, k2, k3) = dGa(k1, k2, k3) + Ga(k1, k6, k4) * xph(k5, k6) * Ga(k5, k2, k3)

        k6 = leg4(k2, k5, k3)
        if(k6 /= 0) dGa(k1, k2, k3) = dGa(k1, k2, k3) + Ga(k1, k6, k5) * xph(k5, k6) * Ga(k2, k5, k3)
    end do  !loop
end do; end do; end do

!----------------------------------------------------



	!update Ga
	dGa = dGa * dL / twopi;    Ga =  Ga + dGa

	if( abs(Vsc) > divergence)exit
	if(abs(Vsdw) > divergence)exit
	if(abs(Vcdw) > divergence)exit


 	!!-----------------analysis
    qsc = 0
    qsdw = pi
    qcdw = pi
    
    call sdw()
    call pdw()
    call cdw()
	write(20,200)Lc, 1/abs(Vsc), 1/abs(Vsdw), 1/abs(Vcdw)
 	
 	
	b11r = 2*Ga(1,2,2)-4*Ga(2,1,2)
	b11s = -(Ga(1,2,2)+1/2*(Ga(1,2,2)-Ga(2,1,2))+1/2*Ga(2,1,2))
	b12r = 2*Ga(3,4,2)-4*Ga(4,3,2)
	b12s = -(Ga(3,4,2)+1/2*(Ga(3,4,2)-Ga(4,3,2))+1/2*Ga(4,3,2))
	f12r = -4*Ga(2,3,2)+2*Ga(3,2,2)
	f12s = 1/2*Ga(2,3,2)+1/2*(Ga(2,3,2)-Ga(3,2,2))-Ga(3,2,2)
	u11r = 4*Ga(2,2,3)+4*Ga(3,3,2)
	u12r = 2*Ga(1,3,4)+2*Ga(2,4,1)+2*Ga(3,1,4)+2*Ga(4,2,1)
	u12s = -Ga(2,4,1)-Ga(3,1,4)

 	write(10, 100)Lc, 1/abs(b11r), 1/abs(b11s), 1/abs(b12r), 1/abs(b12s), 1/abs(f12r), 1/abs(f12s), 1/abs(u11r), 1/abs(u12r), 1/abs(u12s)
 	
 	!!!save data
 	
 	!marki = 2+2*extra; markj = 3+4*extra; markk = 4+6*extra
    
    
    !Gf1 = Ga(markj, markj, markj); Gf2 = Ga(markj, 1, 1)
    !Gb1 = Ga(markj, marki, markj); Gb2 = Ga(markk, markj, marki); Gb3 = Ga(1, markj, 1)
    !Gb4 = Ga(markk, markj, markk); Gb5 = Ga(1, marki, 1)
    !Gu1 = Ga(markk, marki, 1); Gu2 = Ga(marki, markk, 1); Gu3 = Ga(markk, markk, 1)

    !write(10, 100)Lc, Gf1, Gf2, Gb1, Gb2, Gb3, Gb4, Gb5, Gu1, Gu2, Gu3
    

    

end do
100 format(1x, 10e15.7)
200 format(1x, 4e15.7)

contains
	subroutine pdw()
	  implicit none
	  integer :: npatch = nK
	  
      V = 0
      do i1 = 1, npatch;  
      i2 = findpatch(i1, qsc, 1, nK, Kvec);  if(i2 < 0 )cycle;  
         do i4 = 1, npatch; 
         i3 = findpatch(i4, qsc, 1, nK, Kvec);  if(i3 < 0 )cycle;
	        V(i1, i4) = V(i1, i4) + Ga(i1, i2, i3)
      end do;  end do
      V = V / npatch
      
      call check(npatch, V);   V = - V
      call DHEIGEN(npatch, V, eval);   Vsc = eval(ival)
	  return
	end subroutine pdw
	
	subroutine sdw()
      implicit none
      integer :: npatch = nK
      
	  V = 0
      do i1 = 1, npatch; 
      i3 = findpatch(i1, qsdw, -1, nK, Kvec);  if(i3 < 0)cycle;  
         do i4 = 1, npatch;  
         i2 = findpatch(i4, qsdw, -1, nK, Kvec);  if(i2 < 0)cycle; 
            V(i1, i4) = V(i1, i4) + Ga(i1, i2, i3)
            !print*, V
      end do;  end do
      V = V / npatch

      call check(npatch, V)
      call DHEIGEN(npatch, V, eval);  Vsdw = eval(ival)
	  return
	end subroutine sdw

    subroutine cdw()
      implicit none
      integer :: npatch = nK
      
	  V = 0
      do i1 = 1, npatch; 
      i4 = findpatch(i1, qcdw, -1, nK, Kvec);  if(i4 < 0)cycle; 
         do i3 = 1, npatch;  
         i2 = findpatch(i3, qcdw, -1, nK, Kvec);  if(i2 < 0)cycle;  
            V(i1, i3) = V(i1, i3) - Ga(i1, i2, i3) * 2
      end do;  end do

      do i1 = 1, npatch;  
      i3 = findpatch(i1, qcdw, -1, nK, Kvec);  if(i3 < 0)cycle;  
         do i4 = 1, npatch;  
         i2 = findpatch(i4, qcdw, -1, nK, Kvec);  if(i2 < 0)cycle; 
            V(i1, i4) = V(i1, i4) + Ga(i1, i2, i3)
      end do;  end do
    
      V = V / npatch

      call check(npatch, V)
      call DHEIGEN(npatch, V, eval);  Vcdw = eval(ival)

	  return
	end subroutine cdw









end  !end of program





function findpatch(i,Q,info,npatch,kv)
  integer i, info, findpatch, npatch
  real*8, dimension(npatch) :: kv
  real*8 :: twopi = asin(1.d0) * 4
  real*8 Q

  integer k
  real*8 del

  findpatch = -1

  do k = 1, npatch
     if( abs( mod( kv(i) + info * kv(k) - Q, twopi ) ) < 1.e-5 )then
	   findpatch = k; return
	 end if
  end do

  return
end function findpatch

!calculate lindhard function
function lindhard_1d(La, Ek1, vK1, Ek2, vK2, radius)
  implicit none
  real*8 La, Ek1, vK1, Ek2, vK2, radius
  real*8 lindhard_1d
  real*8 :: twopi = asin(1.d0) * 4
  complex*16 z1, z2;  complex*16 :: one = (0, 1)

  !purpose:
  ! integration over momentum as follows
  !    chi = it (dq/2 pi) G(i*La - (Ek1 + vk1 * q)) G(i*La - (Ek2 + vK2 * q) ) + c. c.
  ! since q is measured in units of pi, in the above pi should be set to unity. the 
  ! factor of 1/2 is cancelled by complex conjugation.

  z1 = (one * La - Ek1) / vK1
  z2 = (one * La - Ek2) / vK2

  if( abs(z1 - z2) < 1.d-10 ) then
    lindhard_1d = 2./twopi * real(-2 * radius / (radius * radius - z1 * z1) / ( vK1 * vK1 ))
  else
    lindhard_1d = 2./twopi * real(( log( (radius-z1)/(-radius-z1) ) - log( (radius-z2)/(-radius-z2) ) ) / (z1-z2) / ( vK1 * vK2 ))
  end if
  
  return
end function lindhard_1d


SUBROUTINE DHEIGEN(N,A,EVAL)
  IMPLICIT NONE
  INTEGER N,INFO
  REAL*8 A(N,N),EVAL(N),WORK(3*N)

  IF(N<=0)STOP 'INVALID N @ DHEIGEN'
  IF(N==1)THEN; EVAL(1)=A(1,1); A=1; RETURN; END IF

  CALL DSYEV('V','U',N,A,N,EVAL,WORK,3*N,INFO)
  IF(INFO/=0)STOP 'DHEIGEN FAILED'
  RETURN
END SUBROUTINE DHEIGEN

subroutine check(ndim,V)
  implicit none
  integer ndim
  real*8 V(ndim,ndim)

  integer i,j

  do i = 1, ndim;  do j = 1, i-1
     if( abs( V(i, j) - V(j, i) ) > 1.e-6 ) print*,'hermiticity violated'
  end do;  end do
  
  return
end subroutine check














