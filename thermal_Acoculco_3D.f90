!
!	therm_Ac_3D.f90
!
!	Created by Fernando J. Guerrero-Martinez (fguerrero@igeofisica.unam.mx) on April/18.
!	Finite Volume scheme to solve the steady-state heat transfer equation
!	for the Acoculco caldera, Puebla, Mexico.
!	The heat input for the caldera is based on an estimated Curie isotherm and distributed
!	shallow heat sources associated with recent magmatism in the area.
!====================================================================================================

program acoculco
use omp_lib
implicit none

!Variables definitions

integer i, j, k, ii, jj, ni, nx, ny, nz, nref, max_iter
real x0, xl, y0, yl, z0, zl, dxw, dys, dzb, dxe, dyn, dzt, se, sw, sn, ss, st, sb
real time1, time2, a, b !a and b temporal variables
real tolerance, residual, rad, dpth, d_log, tinter, dev, tamb, tcur, tdike 
character*50 txt, vol, dph, mcase, Filename, frmt 

!Arrays definitions
real, allocatable, dimension(:) :: xc, yc, zc, x, y, z !mesh arrays
real, allocatable, dimension(:,:) :: topo, curie_base, tlogs !Temperatures read from input files 
real, allocatable, dimension(:,:,:) :: t, ap, ae, aw, an, as, at, ab, sp, gamma !Matrix coefficients and variables
integer, allocatable, dimension(:,:,:) :: mark_cells 

!Number of mesh elements
nx=145
ny=167
nz=364

allocate (t(0:nx+1,0:ny+1,0:nz+1), ap(nx,ny,nz), ae(nx,ny,nz), aw(nx,ny,nz))
allocate (an(nx,ny,nz), as(nx,ny,nz), at(nx,ny,nz), ab(nx,ny,nz), sp(nx,ny,nz))
allocate (xc(0:nx+1), yc(0:ny+1), zc(0:nz+1), x(0:nx), y(0:ny), z(0:nz))
allocate (gamma(0:nx+1,0:ny+1,0:nz+1), topo(0:nx+1,0:ny+1), curie_base(0:nx+1,0:ny+1))
allocate (mark_cells(0:nx+1,0:ny+1,0:nz+1), tlogs(43,2))

!set number of threads por parallel computing of residual
call omp_set_num_threads(2)

!Domain definition, UTM coordinates
x0=583340.0
xl=596390.0
y0=2196992.0
yl=2212022.0
z0=-6000.0
zl=3100.0

!maximum iteretions for the solver
max_iter=8000

time1=omp_get_wtime()

!builds the 3D mesh
call Mesh1D (xc, x, x0, xl, nx)
call Mesh1D (yc, y, y0, yl, ny)
call Mesh1D (zc, z, z0, zl, nz)

open(unit=1, file='in_elevAcoculco.txt')
open(unit=2, file='in_CurieDepth.txt')
do i=1, 145
	do j=1, 167
	read(1,*) topo(i,j)
	read(2,*) curie_base(i,j)
	end do
end do
close(1)
close(2)

!extrapolation of topography at boundary nodes
topo(0,:)=topo(1,:); topo(nx+1,:)=topo(nx,:)
topo(:,0)=topo(:,1); topo(:,ny+1)=topo(:,ny)

curie_base(0,:)=curie_base(1,:); curie_base(nx+1,:)=curie_base(nx,:)
curie_base(:,0)=curie_base(:,1); curie_base(:,ny+1)=curie_base(:,ny)

!Set number of depths of the heat sources to evaluate
do ni=1, 6

!=======================================================================
!                    Initial values of variables
t=0.0;aP=0.0;ae=0.0;aw=0.0;an=0.0;as=0.0;at=0.0;ab=0.0;sp=0.0
mark_cells=0

tamb=15.0; tcur=580.0

!=======================================================================
!            Initial temperature field. Geothermal gradient
do i=0, nx+1
	do j=0, ny+1
		do k=0, nz+1
			if(mark_cells(i,j,k) .eq. 1 .or. mark_cells(i,j,k) .eq. 2) cycle
			t(i,j,k)=tamb+(topo(i,j)-zc(k))*0.08
		end do
	end do
end do

!=======================================================================
!	Mark atmospheric control volumes with integer 1 
!	(max elevation 3092masl, min elevation 2261 masl)
do i=0, nx+1
do j=0, ny+1
    do k=0, nz+1
        if(zc(k) .ge. topo(i,j))then
            mark_cells(i,j,k:nz+1)=1
	    t(i,j,k:nz+1)=tamb
            exit
        end if
    end do
end do
end do

!mark control volumes below the Curie isotherm with integer 2
do i=0, nx+1
do j=0, ny+1
    do k=0, nz+1
        if(zc(k) .gt. curie_base(i,j))then
       		mark_cells(i,j,0:k-1)=2
		t(i,j,0:k-1)=tcur			
            exit
        end if
    end do
end do
end do

!=======================================================================
!      Mark 5 cylindrical heat sources at constant temperature 750 C
dpth=-200.0-(100.0*(ni-1)) ! m asl
tdike=750.0
rad=400.0
call select_disc_3D(1, t, mark_cells, 3, x, xc, y, yc, z, zc, nx, ny, nz, rad, 586000.0, 2202840.0, dpth, tdike)
call select_disc_3D(1, t, mark_cells, 3, x, xc, y, yc, z, zc, nx, ny, nz, rad, 587500.0, 2201220.0, dpth, tdike)
call select_disc_3D(1, t, mark_cells, 3, x, xc, y, yc, z, zc, nx, ny, nz, rad, 590000.0, 2200600.0, dpth, tdike)
call select_disc_3D(1, t, mark_cells, 3, x, xc, y, yc, z, zc, nx, ny, nz, rad, 589700.0, 2203000.0, dpth, tdike)
call select_disc_3D(1, t, mark_cells, 3, x, xc, y, yc, z, zc, nx, ny, nz, rad, 589928.0, 2204890.0, dpth, tdike)

!=======================================================================
!       Define thermal conductivity and store values in array gamma
call conductivity_Ac (xc, yc, zc, nx, ny, nz, gamma)

!Additionally, set a conductivity of 0.8 W/(m-K) for the upper 50 m soils
do k=0, nz+1
   do i=0, nx+1
    do j=0, ny+1
        if(mark_cells(i,j,k) .eq. 1 .and. mark_cells(i,j,k-1) .eq. 0)then
        gamma(i,j,k-2:k-1)=0.8 !La superficie	
        end if
    end do
    end do
end do


!Define coefficients for the matrix
	do i=1, nx
		do j=1, ny
			do k=1, nz
			dxw=xc(i)-xc(i-1)
			dys=yc(j)-yc(j-1)
			dzb=zc(k)-zc(k-1)
			
			dxe=xc(i+1)-xc(i)
			dyn=yc(j+1)-yc(j)
			dzt=zc(k+1)-zc(k)
	
			se=(y(j)-y(j-1))*(z(k)-z(k-1))
			sw=(y(j)-y(j-1))*(z(k)-z(k-1))
			sn=(x(i)-x(i-1))*(z(k)-z(k-1))
			ss=(x(i)-x(i-1))*(z(k)-z(k-1))
			st=(x(i)-x(i-1))*(y(j)-y(j-1))
			sb=(x(i)-x(i-1))*(y(j)-y(j-1))

			ae(i,j,k)=0.5*(gamma(i,j,k)+gamma(i+1,j,k))*se/dxe
			aw(i,j,k)=0.5*(gamma(i,j,k)+gamma(i-1,j,k))*sw/dxw
			an(i,j,k)=0.5*(gamma(i,j,k)+gamma(i,j+1,k))*sn/dyn
			as(i,j,k)=0.5*(gamma(i,j,k)+gamma(i,j-1,k))*ss/dys
			at(i,j,k)=0.5*(gamma(i,j,k)+gamma(i,j,k+1))*st/dzt
			ab(i,j,k)=0.5*(gamma(i,j,k)+gamma(i,j,k-1))*sb/dzb
			ap(i,j,k)=ae(i,j,k)+aw(i,j,k)+an(i,j,k)+as(i,j,k)+at(i,j,k)+ab(i,j,k)
			sp(i,j,k)=0.0
		
		!The surface, Curie isotherm, and heat sources are kept at a constant temperature
		!the following conditions take this into account and correct coefficients accordingly
		if (mark_cells(i,j,k) .eq. 1 .or. mark_cells(i,j,k) .eq. 2 .or. mark_cells(i,j,k) .eq. 3) then
		ae(i,j,k)=0.0; aw(i,j,k)=0.0; an(i,j,k)=0.0; as(i,j,k)=0.0; at(i,j,k)=0.0; ab(i,j,k)=0.0; ap(i,j,k)=1.0
			if(mark_cells(i,j,k) .eq. 1) sp(i,j,k)=tamb
			if(mark_cells(i,j,k) .eq. 2) sp(i,j,k)=tcur
			if(mark_cells(i,j,k) .eq. 3) sp(i,j,k)=tdike
		else
		
			if (mark_cells(i+1,j,k) .eq. 1 .or. mark_cells(i+1,j,k) .eq. 2 .or. mark_cells(i+1,j,k) .eq. 3) then
				ap(i,j,k)=ap(i,j,k)+ae(i,j,k)
				if (mark_cells(i+1,j,k) .eq. 1) sp(i,j,k)=sp(i,j,k)+2.0*ae(i,j,k)*tamb
			        if (mark_cells(i+1,j,k) .eq. 2) sp(i,j,k)=sp(i,j,k)+2.0*ae(i,j,k)*tcur
				if (mark_cells(i+1,j,k) .eq. 3) sp(i,j,k)=sp(i,j,k)+2.0*ae(i,j,k)*tdike	 
				ae(i,j,k)=0.0
			end if
			
			if (mark_cells(i-1,j,k) .eq. 1 .or. mark_cells(i-1,j,k) .eq. 2 .or. mark_cells(i-1,j,k) .eq. 3) then
				ap(i,j,k)=ap(i,j,k)+aw(i,j,k)
				if (mark_cells(i-1,j,k) .eq. 1) sp(i,j,k)=sp(i,j,k)+2.0*aw(i,j,k)*tamb
				if (mark_cells(i-1,j,k) .eq. 2) sp(i,j,k)=sp(i,j,k)+2.0*aw(i,j,k)*tcur
				if (mark_cells(i-1,j,k) .eq. 3) sp(i,j,k)=sp(i,j,k)+2.0*aw(i,j,k)*tdike 
				aw(i,j,k)=0.0
			end if

			if (mark_cells(i,j+1,k) .eq. 1 .or. mark_cells(i,j+1,k) .eq. 2 .or. mark_cells(i,j+1,k) .eq. 3) then
				ap(i,j,k)=ap(i,j,k)+an(i,j,k)
				if (mark_cells(i,j+1,k) .eq. 1) sp(i,j,k)=sp(i,j,k)+2.0*an(i,j,k)*tamb
				if (mark_cells(i,j+1,k) .eq. 2) sp(i,j,k)=sp(i,j,k)+2.0*an(i,j,k)*tcur
				if (mark_cells(i,j+1,k) .eq. 3) sp(i,j,k)=sp(i,j,k)+2.0*an(i,j,k)*tdike
				an(i,j,k)=0.0
			end if

			if (mark_cells(i,j-1,k) .eq. 1 .or. mark_cells(i,j-1,k) .eq. 2 .or. mark_cells(i,j-1,k) .eq. 3) then
				ap(i,j,k)=ap(i,j,k)+as(i,j,k)
				if (mark_cells(i,j-1,k) .eq. 1) sp(i,j,k)=sp(i,j,k)+2.0*as(i,j,k)*tamb
				if (mark_cells(i,j-1,k) .eq. 2) sp(i,j,k)=sp(i,j,k)+2.0*as(i,j,k)*tcur
				if (mark_cells(i,j-1,k) .eq. 3) sp(i,j,k)=sp(i,j,k)+2.0*as(i,j,k)*tdike 
				as(i,j,k)=0.0
			end if
			
			if (mark_cells(i,j,k+1) .eq. 1 .or. mark_cells(i,j,k+1) .eq. 3) then
				ap(i,j,k)=ap(i,j,k)+at(i,j,k)
				if (mark_cells(i,j,k+1) .eq. 1) sp(i,j,k)=sp(i,j,k)+2.0*at(i,j,k)*tamb 
				if (mark_cells(i,j,k+1) .eq. 3) sp(i,j,k)=sp(i,j,k)+2.0*at(i,j,k)*tdike 
				at(i,j,k)=0.0
			end if

			if (mark_cells(i,j,k-1) .eq. 2 .or. mark_cells(i,j,k-1) .eq. 3) then
				ap(i,j,k)=ap(i,j,k)+ab(i,j,k)
				if (mark_cells(i,j,k-1) .eq. 2) sp(i,j,k)=sp(i,j,k)+2.0*ab(i,j,k)*tcur
				if (mark_cells(i,j,k-1) .eq. 3) sp(i,j,k)=sp(i,j,k)+2.0*ab(i,j,k)*tdike 
				ab(i,j,k)=0.0
			end if
		end if

			end do
		end do
	end do

!Boundary conditions

!East
	ap(nx,:,:)=ap(nx,:,:)-ae(nx,:,:)
	ae(nx,:,:)=0.0

!West
	ap(1,:,:)=ap(1,:,:)-aw(1,:,:)
	aw(1,:,:)=0.0

!North
	ap(:,ny,:)=ap(:,ny,:)-an(:,ny,:)
	an(:,ny,:)=0.0

!South
	ap(:,1,:)=ap(:,1,:)-as(:,1,:)
	as(:,1,:)=0.0

!Top
	sp(:,:,nz)=sp(:,:,nz)+at(:,:,nz)*tamb
	at(:,:,nz)=0.0

!Bottom
	sp(:,:,1)=sp(:,:,1)+ab(:,:,1)*tcur
	ab(:,:,1)=0.0

	tolerance=1.5E-2 !convergence criterion
	call Gauss_Seidel3D(t,nx,ny,nz,aP,ae,aw,an,as,at,ab,sp,nx,ny,nz,max_iter,tolerance,residual)	

!=======================================================================
!			OUTPUT
!		 3D temperature field
frmt='(F8.1, A, F9.1, A, F7.1, A, F6.1)'
write(mcase,'(i6)') int(abs(dpth))
mcase=ADJUSTL(mcase)
Filename='temp3D_'//mcase(1:len_trim(mcase))//".txt"
open(1,file=Filename(1:len_trim(Filename)))

write (1,*) 'X',', ', 'Y',', ', 'Z',', ', 'T'
    do i=1, nx
        do j=1, ny
            do k=120, nz
    if(mod(float(k),2.0) .eq. 0.0) then
        if(mark_cells(i,j,k) .eq. 1) then
        write(1,frmt) xc(i),', ', yc(j),', ', zc(k),', ',  2.0
        else
        write(1,frmt) xc(i),', ', yc(j),', ', zc(k),', ',  t(i,j,k)
        end if
    end if
            end do
        end do
    end do
close(1)

! 		2D section of Temperature
nref=73
	Filename='temp2D_'//mcase(1:len_trim(mcase))//".txt"
	open(1,file=Filename(1:len_trim(Filename)))
do j=1, ny
	do k=1, nz
		if(mark_cells(nref,j,k) .eq. 1) then
			write(1,*) yc(j), zc(k), 2.0
		else
			write(1,*) yc(j), zc(k), t(nref,j,k)
		end if
	end do
	write(1,*) ' '
end do
close(1)


!=======================================================================  
!Write thermal profile of EAC-1 and its overall deviation from the temperature log
frmt='(F7.1, F6.1)'
! EAC-1
! get the mesh nodes of the well
call get_mesh_node(x, nx, 589700.0, ii, 1)
call get_mesh_node(y, ny, 2203000.0, jj, 1)

Filename='eac1_'//mcase(1:len_trim(mcase))//".txt"

! the thermal profile is printed up to 3000 m depth
open(1,file=Filename(1:len_trim(Filename)))
do k=nz+1, 0, -1
	if(topo(ii,jj)-zc(k) .lt. 0.0) cycle
	if(topo(ii,jj)-zc(k) .gt. 3000.0) exit
	write(1,frmt) -(topo(ii,jj)-zc(k)), t(ii,jj,k) 
end do
close(1)

!============== Overall deviation from measured temperatures
!read temperature logs of EAC-1 and EAC-2
!in_tLogs: EAC-1 from index 1 to 18; EAC-2 from 19 owards
open(unit=1, file='in_tLogs.txt')
do i=1, 42
	read(1,*) tlogs(i,1), tlogs(i,2)
end do
close(1)

dev=0.0
do i=1, 18
d_log=tlogs(i,2)

k=nz
do while (topo(ii,jj)-zc(k) .lt. d_log)
k=k-1
end do
a=topo(ii,jj)-zc(k+1); b=topo(ii,jj)-zc(k)

tinter=t(ii,jj,k+1)+((d_log-a)/(b-a))*(t(ii,jj,k)-t(ii,jj,k+1))

dev=dev+abs(tinter-tlogs(i,1))
end do
!write overall deviation of well EAC-1
Filename='deviation_eac1_'//mcase(1:len_trim(mcase))//".txt"
open(1,file=Filename(1:len_trim(Filename)))
write(1,*) dev/17.0
close(1)

! EAC-2. Only temperature profile is printed
call get_mesh_node(x, nx, 590172.0, ii, 1)
call get_mesh_node(y, ny, 2203254.0, jj, 1)

Filename='eac2_'//mcase(1:len_trim(mcase))//".txt"
open(1,file=Filename(1:len_trim(Filename)))
do k=nz+1, 0, -1
	if(topo(ii,jj)-zc(k) .lt. 0.0) cycle
	if(topo(ii,jj)-zc(k) .gt. 3000.0) exit
	write(1,frmt) -(topo(ii,jj)-zc(k)), t(ii,jj,k) 
end do
close(1)



end do


time2=omp_get_wtime()
write(*,*)  'cpu time=', time2-time1

end program

!=======================================================================
!Set the thermal conducivity in the 3D domain
Subroutine conductivity_Ac (xc, yc, zc, nx, ny, nz, gamma)
implicit none
integer i, j, k, nx, ny, nz
real*8 bound
real xc(0:nx+1), yc(0:ny+1), zc(0:nz+1), gamma(0:nx+1,0:ny+1,0:nz+1), rho(nz), cesp(nz)

!write(*,*)'---> Conductivities'
do k=0, nz+1
    do i=0, nx+1
    do j=0, ny+1

if(yc(j) .le. -0.5381*xc(i)+2520592.0)then
        if(zc(k) .le. 1240.0)then!Below this depth there is either granite or limestone and marl
            bound=1091.0-1.9E-4*((xc(i)-589700)**2+(yc(j)-2203000)**2)
            if(zc(k) .le. bound)then
		gamma(i,j,k)=2.20 !granite
            else
                gamma(i,j,k)=2.80 !Limestone
            end if
        end if

        if(zc(k) .ge. 1240.0 .and. zc(k) .lt. 1362.0)then
		if(yc(j) .ge. -0.5381*xc(i)+2520592.0-450.0 .and. (xc(i) .lt. 589700.0+450.0 .and. xc(i) .gt. 589700.0-450.0 ))then
                   gamma(i,j,k)=2.90 !Aplite
		else
		   gamma(i,j,k)=2.60 !Limestone	
		end if
        end if

        if(zc(k) .ge. 1362.0 .and. zc(k) .lt. 1575.0)then
		if(yc(j) .ge. -0.5381*xc(i)+2520592.0-450.0 .and. (xc(i) .lt. 589700.0+450.0 .and. xc(i) .gt. 589700.0-450.0 ))then
                   gamma(i,j,k)=1.60 !Skarn
		else
		   gamma(i,j,k)=2.60 !Limestone	
		end if
        end if

        if(zc(k) .ge. 1575.0 .and. zc(k) .lt. 1670.0)then
        gamma(i,j,k)=2.80 !marl
        end if

        if(zc(k) .ge. 1670.0 .and. zc(k) .lt. 2040.0)then
		if(yc(j) .ge. -0.5381*xc(i)+2520592.0-450.0 .and. (xc(i) .lt. 589700.0+450.0 .and. xc(i) .gt. 589700.0-450.0 ))then
                   gamma(i,j,k)=1.60 !Skarn
		else
		   gamma(i,j,k)=2.60 !Limestone	
		end if
        end if

        if(zc(k) .ge. 2040.0 .and. zc(k) .lt. 2581.0)then
        gamma(i,j,k)=2.08 !Riodacite
        end if

        if(zc(k) .ge. 2581.0 .and. zc(k) .lt. 2646.0)then
        gamma(i,j,k)=1.65 !Ignimbrite
        end if

        if(zc(k) .ge. 2646.0 .and. zc(k) .lt. 2736.0)then
        gamma(i,j,k)=1.87 !Dacite
        end if

        if(zc(k) .ge. 2736.0)then
        gamma(i,j,k)=1.40 !Tuffs
        end if
else !Lithology of well EAC2
        if(zc(k) .lt. 2319.0)then! Below this depth there are "Hornfels", "Granite" and "Limestone"
            bound=1204.0-1.9E-4*((xc(i)-590172.0)**2+(yc(j)-2203254.0)**2)
            if(zc(k) .le. bound)then
		gamma(i,j,k)=2.20!Granite
            else
		if(yc(j) .le. -0.5381*xc(i)+2520592.0+450.0 .and. (xc(i) .lt. 590172.0+450.0 .and. xc(i) .gt. 590172.0-450.0 ))then
                   gamma(i,j,k)=2.50 !Hornfels
		else
		   gamma(i,j,k)=2.60 !Limestone	
		end if
            end if
        end if

        if(zc(k) .ge. 2319.0 .and. zc(k) .lt. 2416.0)then
                gamma(i,j,k)=2.60 !Limestone
        end if

        if(zc(k) .ge. 2416.0 .and. zc(k) .lt. 2479.0)then
        gamma(i,j,k)=1.72 !Andesite
        end if

        if(zc(k) .ge. 2479.0 .and. zc(k) .lt. 2606.0)then
        gamma(i,j,k)=1.87 !Dacite
        end if

        if(zc(k) .ge. 2606.0 .and. zc(k) .lt. 2639.0)then
        gamma(i,j,k)=1.72 !Andesite
        end if

        if(zc(k) .ge. 2639.0)then
        gamma(i,j,k)=1.65 !Tuff
        end if


end if

    end do
    end do
end do

end subroutine

!=======================================================================
!Transforms distance to a number of node 
!control: 1) returns node of the control volume; 
!	  2) returns number of face of the control volume
!=======================================================================
subroutine get_mesh_node(x, nx, xi, nvc, control)
implicit none
integer nx, nvc, control
real  x0, xl, xi
real  x(0:nx)

x0=x(0); xl=x(nx)

	if(xi<x0 .or. xi>xl)then 
		write(*,*) 'get_mesh_node: distance out of the domain', xi, x0, xl
		stop
	end if

select case (control)
!case 1, returns node
	case (1)
		nvc=0
		do while(x(nvc)<xi)
			nvc=nvc+1
		end do

		if(xi .eq. x(nx) )then
			nvc=nvc+1
		end if
!case 2, returns face
	case (2)
		nvc=0
		do while(x(nvc)<xi)
			nvc=nvc+1
		end do
		
		if((xi .ne. x(0)).and. (xi-x(nvc-1)<x(nvc)-xi))then
			nvc=nvc-1
		end if
end select
end subroutine

!=======================================================================
!Asigns a temperature value to a disc of a specified radius on on a plane x, y, or z 
!Control: 1, 2, and 3 for planes z, x, y respectively.
!Depth and radius are specified
!marker is an iteger to identify the location of the disc
!=======================================================================

subroutine select_disc_3D(control, t, mark_cells, marker, x, xc, y, yc, z, zc, nx, ny, nz, rad, xi, xj, depth, temp)

implicit none
integer  i, j, k, nx, ny, nz, control, n_depth, ncen1, ncen2, marker
real  x0, xl, y0, yl, z0, zl, xi, xj, rad, temp, depth
real  xc(0:nx+1), yc(0:ny+1), zc(0:nz+1), x(0:nx), y(0:ny), z(0:nz), t(0:nx+1,0:ny+1,0:nz+1)
integer mark_cells(0:nx+1,0:ny+1,0:nz+1)

x0=xc(0); xl=xc(nx+1); y0=yc(0); yl=yc(ny+1); z0=zc(0); zl=zc(nz+1)

SELECT CASE (control)
	CASE (1) !plane z constant
		call get_mesh_node(x, nx, xi, ncen1, 2)
		call get_mesh_node(y, ny, xj, ncen2, 2)
		call get_mesh_node(z, nz, depth, n_depth, 1)
			do i=1, nx
				do j=1, ny
					if(sqrt((xc(i)-x(ncen1))**2+(yc(j)-y(ncen2))**2)<rad)then
						t(i,j,n_depth-1:n_depth)=temp
						mark_cells(i,j,n_depth-1:n_depth)=marker
					end if
				end do
			end do
	CASE (2)!plane x constant
		call get_mesh_node(y, ny, xi, ncen1, 2)
		call get_mesh_node(z, nz, xj, ncen2, 2)
		call get_mesh_node(x, nx, depth, n_depth, 1)
			do j=1, ny
				do k=1, nz
					if(sqrt((yc(j)-y(ncen1))**2+(zc(k)-z(ncen2))**2)<rad)then
						t(n_depth,j,k)=temp
						mark_cells(n_depth,j,k)=marker
					end if
				end do
			end do
	CASE (3)!plane y constant
		call get_mesh_node(x, nx, xi, ncen1, 2)
		call get_mesh_node(z, nz, xj, ncen2, 2)
		call get_mesh_node(y, ny, depth, n_depth, 1)

			do i=1, nx
				do k=1, nz
					if(sqrt((xc(i)-x(ncen1))**2+(zc(k)-z(ncen2))**2)<rad)then
						t(i,n_depth,k)=temp
						mark_cells(i,n_depth,k)=marker
					end if
				end do
			end do
end select
end subroutine

!=======================================================================
!Builds a uniform mesh. xc stores node and x sotres face of control volumes
Subroutine Mesh1D(xc,x,x0,xl,nx)
integer i,j,nx
real*4 x0,xl,dx
real*4 x(0:nx),xc(0:nx+1)
dx=(1.0)/dfloat(nx)

do i=0,nx
	x(i)=dfloat(i)*dx
	x(i)=x0+(xl-x0)*x(i)
end do


xc(0)=x(0); xc(nx+1)=x(nx)
do i=1,nx
xc(i)=(x(i)+x(i-1))*0.5
end do

End Subroutine

!=======================================================================
!Solver
subroutine Gauss_Seidel3D(phi,ei,ej,ek,aP,aE,aW,aN,aS,aF,aB,sP,nx,ny,nz,max_iter,tolerance,residual)

implicit none

integer ei,ej,ek,i,j,k,nx,ny,nz,count_iter,max_iter
real*4 phi(0:ei+1,0:ej+1,0:ek+1),aP(nx,ny,nz),aE(nx,ny,nz),aW(nx,ny,nz),aN(nx,ny,nz), &
aS(nx,ny,nz),aF(nx,ny,nz),aB(nx,ny,nz),sP(nx,ny,nz)
real*4 residual,tolerance

interface
real function calcResidual(phi,ei,ej,ek,aP,aE,aW,aN,aS,aF,aB,sP,nx,ny,nz)
implicit none
integer ei,ej,ek,nx,ny,nz
real phi(0:ei+1,0:ej+1,0:ek+1),aP(nx,ny,nz),aE(nx,ny,nz),aW(nx,ny,nz),aN(nx,ny,nz), &
aS(nx,ny,nz),aF(nx,ny,nz),aB(nx,ny,nz),sP(nx,ny,nz)
end function
end interface


count_iter=0;  residual=1.0

do while((count_iter <= max_iter).and.(residual > tolerance))
do i=1,ei
do j=1,ej
do k=1,ek
phi(i,j,k)=(aE(i,j,k) * phi(i+1,j,k) + aW(i,j,k) * phi(i-1,j,k)+aN(i,j,k) * phi(i,j+1,k) + &
aS(i,j,k) * phi(i,j-1,k)+aB(i,j,k) * phi(i,j,k-1)+aF(i,j,k) * phi(i,j,k+1)+sp(i,j,k))/aP(i,j,k)
end do
end do
end do
residual = calcResidual(phi,ei,ej,ek,aP,aE,aW,aN,aS,aF,aB,sP,nx,ny,nz)
write(*,*) 'Residual', count_iter, residual
count_iter=count_iter+1
end do

end subroutine

!=======================================================================
!Residual
real function calcResidual(phi,ei,ej,ek,aP,aE,aW,aN,aS,aF,aB,sP,nx,ny,nz)
implicit none
real residual,NINV
integer bi,ei,bj,ej,bk,ek,i,j,k,nx,ny,nz
real phi(0:ei+1,0:ej+1,0:ek+1),aP(nx,ny,nz),aE(nx,ny,nz),aW(nx,ny,nz),aN(nx,ny,nz), &
aS(nx,ny,nz),aF(nx,ny,nz),aB(nx,ny,nz),sP(nx,ny,nz)
real, allocatable :: acum(:,:,:)
allocate (acum(ei,ej,ek))

bi=1; bj=1; bk=1.0
acum=0
NINV = 1.0 / dfloat(ei*ej*ek)

do i=bi,ei
do j=bj,ej
	!$omp parallel
	!$omp do schedule(static)
do k=bk,ek

acum(i,j,k) = (aE(i,j,k) * phi(i+1,j,k) +aW(i,j,k) * phi(i-1,j,k) +&
	aN(i,j,k) * phi(i,j+1,k) +aS(i,j,k) * phi(i,j-1,k)+aF(i,j,k) * phi(i,j,k+1) &
	+aB(i,j,k) * phi(i,j,k-1)+sP(i,j,k))-aP(i,j,k) * phi(i,j,k)

end do
	!$omp end parallel
end do
end do

residual = sqrt( NINV * sum(acum * acum) )
calcResidual=residual
end function

!=======================================================================

Subroutine computing_time
real*4 cputime, min, hou, sec

call cpu_time(cputime)

if(cputime .lt. 60.0)then
write(*,*)'computing time:', cputime, 'seconds'
end if

if((cputime .ge. 60.0) .and. (cputime .lt. 3600.0))then
min=int(cputime/60.0)
sec=(cputime/60.0-int(cputime/60.0))*60.0
write(*,*)'computing time:', min, 'minutes', sec, 'seconds'
end if

if(cputime .ge. 3600.0)then
hou=int(cputime/3600.0)
min=int((cputime/3600.0-int(cputime/3600.0))*60.0)
sec=(((cputime/3600.0-int(cputime/3600.0))*60.0)-int((cputime/3600.0-int(cputime/3600.0))*60.0))*60.0
write(*,*)'computing time:', hou, 'hours,', min, 'minutes,', sec, 'seconds'
end if

end Subroutine
