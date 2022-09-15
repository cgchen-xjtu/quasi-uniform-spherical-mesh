# quasi-uniform-spherical-mesh
Fortran codes to build quasi-uniform spherical mesh


GnomonicCubedSphereGRID.f90
Fortran codes to build gnomonic cubed sphere 

Spcify the resolution by changing parameter N

NX=(N-1)/2, NY=(N-1)/2 are no. of computational cells used in multi-moment scheme

Function list:

subroutine pprop2sp(lambda,theta,x,y,k)

find the location in global lon-lat coordinates

INPUT: 
(x,y) the location in local coordinates on patch k

subroutine pprosp2p(x,y,lambda,theta,k)

find the location in local coordinates on patch k

INPUT: 
(lambda,theta) the location in global lon-lat coordinates
k the index of the targeted patch

subroutine contravprosp2p(contrav1,contrav2,sv1,sv2,k,lambda,theta)

calculate the contravariant velocity components in local coordinates on patch k from the velocty vector in global lon-lat coordinates

INPUT: 
(sv1,sv2) the velocty vector in global lon-lat coordinates
(lambda,theta) the location in global lon-lat coordinates
k the index of the targeted patch

subroutine contravprop2sp(sv1,sv2,contrav1,contrav2,k,lambda,theta)

calculate the velocty vector in global lon-lat coordinates from the contravariant velocity components in local coordinates on patch k

INPUT: (contrav1,contrav2) the contravariant velocity components in local coordinates on patch k
(lambda,theta) the location in global lon-lat coordinates
k the index of the targeted patch

subroutine matrixA(ma,k,lambda,theta)

subroutine matrixIA(ima,k,lambda,theta)

calculate the base vector and its inverse on patch k

INPUT: 
(lambda,theta) the location in global lon-lat coordinates
k the index of the targeted patch

subroutine matrixG(mg,xx,yy)

subroutine matrixIG(img,xx,yy)

calculate the contravariant and covariant metric tensors

INPUT: 
(x,y) the location in local coordinates on patch k

subroutine computejab(jab,xx,yy)

calculate the Jocabian

INPUT: 
(x,y) the location in local coordinates on patch k

subroutine ghostlocation()

find the location of ghost cells

icosahedral-hexagonalGRID.f90
Fortran codes to build icosahedral-hexagonal grid

specfy the resolution by changing parameter P

No. of vertices : 20*P*P
No. of lines : 30*P*P
No. of cells : 10*P*P+2

the output file is P.dat

the output file include

			DO I=1,20*P*P 
				WRITE(10,*)NODE(I,:)*RADIUS   ! the location of the vertices
			ENDDO
      
			DO I=1,10*P*P+2
				WRITE(10,*)CELL(I,:)*RADIUS   ! the location of the cell barycenter
			ENDDO
      
			DO I=1,30*P*P 
				WRITE(10,*)NODE_ONLINE(I,:)*RADIUS   ! the location of the midpoint of line 
			ENDDO

			DO I=1,10*P*P+2
				WRITE(10,*)CELL_NODES(I,:)    ! the indices of vertices belonging to a cell
			ENDDO
      
			DO I=1,10*P*P+2
				WRITE(10,*)CELL_LINES(I,:)    ! the indices of lines belonging to a cell
			ENDDO

			DO I=1,20*P*P
				WRITE(10,*)NODE_CELLS(I,:)    ! the indices of cells sharing a vertex
			ENDDO

			DO I=1,20*P*P
				WRITE(10,*)NODE_NODES(I,:)    ! the indices of lines belonging to a cell
			ENDDO

			DO I=1,30*P*P
				WRITE(10,*)LINE_NODES(I,:)    ! the indices of lines belonging to a cell
			ENDDO
      
			DO I=1,30*P*P
				WRITE(10,*)LINE_CELLS(I,:)    ! the indices of 2 cells sharing a line
			ENDDO

