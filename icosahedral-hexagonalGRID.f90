!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	PROGRAM MAIN
		IMPLICIT NONE
		INTEGER I,J,K,P,M 
		CHARACTER(20) CH
		PARAMETER(P=10)
		DOUBLE PRECISION,DIMENSION(1:20*P*P,1:3)::NODE
		DOUBLE PRECISION,DIMENSION(1:10*P*P+2,1:3)::CELL
		DOUBLE PRECISION,DIMENSION(1:30*P*P,1:3)::NODE_ONLINE

		DOUBLE PRECISION,DIMENSION(1:20*P*P,1:2)::NODES
		DOUBLE PRECISION,DIMENSION(1:10*P*P+2,1:2)::CELLS
		DOUBLE PRECISION,DIMENSION(1:30*P*P,1:2)::NODES_ONLINE

		INTEGER ,DIMENSION(1:10*P*P+2,1:6)::CELL_NODES,CELL_LINES
		INTEGER ,DIMENSION(1:20*P*P,1:3)::NODE_CELLS,NODE_CELLS_OLD,NODE_LINES,NODE_LINES_OLD,NODE_NODES
		INTEGER ,DIMENSION(1:30*P*P,1:2)::LINE_NODES,LINE_CELLS

		DOUBLE PRECISION TEST,BETA,RADIUS

		DOUBLE PRECISION,DIMENSION(1:3)::PCELL
!--------------------------------------------
		BETA=DACOS(1.D0/DSQRT(5.D0))
		RADIUS=6371220.D0

		CALL MESH(P,NODE,CELL,NODE_CELLS_OLD,CELL_NODES)		!OK	
		CALL MESH_LINES(P,NODE_CELLS_OLD,CELL_NODES,LINE_NODES,LINE_CELLS,CELL_LINES,NODE_LINES_OLD)
		CALL CHANGE_NODECELLS(P,NODE,CELL,LINE_NODES,NODE_CELLS_OLD,NODE_LINES_OLD,CELL_LINES,NODE_CELLS,NODE_LINES,NODE_NODES)
		CALL MESH_CENTER_POINT(P,NODE,LINE_NODES,NODE_ONLINE)

!		DO I=1,20*P*P 
!			NODES(I,1)=DATAN2(NODE(I,2),NODE(I,1))
!			NODES(I,2)=DASIN(NODE(I,3))
!		ENDDO
!		
!		DO I=1,30*P*P 
!			NODES_ONLINE(I,1)=DATAN2(NODE_ONLINE(I,2),NODE_ONLINE(I,1))
!			NODES_ONLINE(I,2)=DASIN(NODE_ONLINE(I,3))
!		ENDDO 
!
!		DO I=1,10*P*P+2
!			TEST=DABS(DABS(CELL(I,3))-1.D0)
!			IF(TEST.GE.0.1D0**10)THEN
!				CELLS(I,1)=DATAN2(CELL(I,2),CELL(I,1))
!				CELLS(I,2)=DASIN(CELL(I,3))
!			ELSE
!				CELLS(I,1)=0.D0
!				CELLS(I,2)=DASIN(CELL(I,3))-DSIGN(1.D0,DASIN(CELL(I,3)))*BETA/DBLE(20*P)
!
!				CELL(I,1)=DCOS(CELLS(I,1))*DCOS(CELLS(I,2))
!				CELL(I,2)=DSIN(CELLS(I,1))*DCOS(CELLS(I,2))
!				CELL(I,3)=DSIN(CELLS(I,2))
!			ENDIF
!		ENDDO
		
		DO I=1,10*P*P+2
			IF(CELL_NODES(I,1)==CELL_NODES(I,6))THEN
				M=5
			ELSE
				M=6
			ENDIF
			IF(M==6)THEN
				CALL CHANGE_CELL(NODE(CELL_NODES(I,1),:),NODE(CELL_NODES(I,3),:),		&
				&		NODE(CELL_NODES(I,4),:),NODE(CELL_NODES(I,6),:),PCELL)
				CELL(I,:)=PCELL(:)
			ENDIF
		ENDDO		

		OPEN(1,FILE='CELL_NODES_NUMBER.DAT',STATUS='UNKNOWN')		
			DO I=1,10*P*P+2		
				WRITE(1,*)CELL(I,:)
			ENDDO
		CLOSE(1)

		OPEN(6,FILE='CELLNODES_Final.DAT',STATUS='UNKNOWN')		
			WRITE(6,*)'TITLE=OUTPUT ICOSAHEDRAL GRID'
			WRITE(6,*)'VARIABLES="X","Y","Z"'
			WRITE(6,*)'ZONE N=',20*P*P,',	E=',6*(12+30*(P-1)+(P-1)*(P-2)*10),',	ZONETYPE=FELINESEG',',DATAPACKING=POINT'		

			DO I=1,20*P*P
				WRITE(6,*)NODE(I,:)
			ENDDO
			DO K=1,10*P*P+2		
				DO J=1,5
					WRITE(6,*)CELL_NODES(K,J),CELL_NODES(K,J+1)
				ENDDO
				WRITE(6,*)CELL_NODES(K,6),CELL_NODES(K,1)
			ENDDO
		CLOSE(6)

		OPEN(2,FILE='TRIANGLE.DAT',STATUS='UNKNOWN')
		WRITE(2,*)'TITLE   =OUTPUT OF AMR ADVECTION SCHEME'

		WRITE(2,*)'VARIABLES = "X" "Y" "Z" "FINE" '

		WRITE(2,*)'ZONE N=',20*P*P+12+30*(P-1)+(P-1)*(P-2)*10,', E=',6*(12+30*(P-1)+(P-1)*(P-2)*10)-12,', ZONETYPE=FETRIANGLE', ', DATAPACKING=POINT'

		DO I=1,20*P*P
			WRITE(2,*)NODE(I,:),1.D0
		ENDDO
		DO I=1,10*P*P+2
			WRITE(2,*)CELL(I,:),1.D0
		ENDDO
		DO I=1,12+30*(P-1)+(P-1)*(P-2)*10
			WRITE(2,*)CELL_NODES(I,1),CELL_NODES(I,2),I+20*P*P
			WRITE(2,*)CELL_NODES(I,2),CELL_NODES(I,3),I+20*P*P
			WRITE(2,*)CELL_NODES(I,3),CELL_NODES(I,4),I+20*P*P
			WRITE(2,*)CELL_NODES(I,4),CELL_NODES(I,5),I+20*P*P
			WRITE(2,*)CELL_NODES(I,5),CELL_NODES(I,6),I+20*P*P
			IF(CELL_NODES(I,1).NE.CELL_NODES(I,6))THEN
				WRITE(2,*)CELL_NODES(I,6),CELL_NODES(I,1),I+20*P*P
			ENDIF
		ENDDO
		CLOSE(2)
!
!		OPEN(2,FILE='CELL_NODES.DAT',STATUS='UNKNOWN')		
!			DO K=1,10*P*P+2		
!				WRITE(2,*)CELL_NODES(K,:)
!			ENDDO
!		CLOSE(2)
!
!		OPEN(3,FILE='LINE_NODES.DAT',STATUS='UNKNOWN')		
!			WRITE(3,*)'TITLE=OUTPUT ICOSAHEDRAL GRID'
!			WRITE(3,*)'VARIABLES="X","Y","Z"'
!			WRITE(3,*)'ZONE N=',20*P*P,',	E=',30*P*P,',	ZONETYPE=FELINESEG',',DATAPACKING=POINT'		
!
!			DO I=1,20*P*P
!				WRITE(3,*)NODE(I,1),NODE(I,2),NODE(I,3)
!			ENDDO
!			DO I=1,30*P*P				
!				WRITE(3,*)LINE_NODES(I,:)
!			ENDDO
!		CLOSE(3)
		
		CALL INT_CHAR(P,CH)

		OPEN(10,FILE=CH,STATUS='UNKNOWN')
			DO I=1,20*P*P 
				WRITE(10,*)NODE(I,:)*RADIUS
			ENDDO
			DO I=1,10*P*P+2
				WRITE(10,*)CELL(I,:)*RADIUS
			ENDDO
			DO I=1,30*P*P 
				WRITE(10,*)NODE_ONLINE(I,:)*RADIUS
			ENDDO

			DO I=1,10*P*P+2
				WRITE(10,*)CELL_NODES(I,:)
			ENDDO
			DO I=1,10*P*P+2
				WRITE(10,*)CELL_LINES(I,:)
			ENDDO

			DO I=1,20*P*P
				WRITE(10,*)NODE_CELLS(I,:)
			ENDDO
!			DO I=1,20*P*P
!				WRITE(10,*)NODE_LINES(I,:)
!			ENDDO
			DO I=1,20*P*P
				WRITE(10,*)NODE_NODES(I,:)
			ENDDO

			DO I=1,30*P*P
				WRITE(10,*)LINE_NODES(I,:)
			ENDDO
			DO I=1,30*P*P
				WRITE(10,*)LINE_CELLS(I,:)
			ENDDO
		CLOSE(10)
	ENDPROGRAM MAIN
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
SUBROUTINE CHANGE_CELL(P1,P3,P4,P6,PCELL)
	IMPLICIT NONE
	DOUBLE PRECISION,DIMENSION(1:3)::P1,P3,P4,P6,PC,PCELL
	DOUBLE PRECISION,DIMENSION(1:3)::P1_NEW,P3_NEW,P4_NEW,P6_NEW
	DOUBLE PRECISION,DIMENSION(1:3)::XX,YY,ZZ
	DOUBLE PRECISION R,R0

	R=DSQRT(P1(1)**2+P1(2)**2+P1(3)**2)

	CALL FIND_LINE_DIRECTION(P6(:)/R,P3(:)/R,YY)	
	CALL FIND_LINE_NDIRECTION(P6(:)/R,YY,ZZ)
	XX(:)=P6(:)/R

	CALL FIND_NEW_VECTOR(XX,YY,ZZ,P1,P1_NEW)
	CALL FIND_NEW_VECTOR(XX,YY,ZZ,P4,P4_NEW)


	PC(1)=P1_NEW(1)-P1_NEW(3)*(P4_NEW(1)-P1_NEW(1))/(P4_NEW(3)-P1_NEW(3))

	PC(2)=P1_NEW(2)-P1_NEW(3)*(P4_NEW(2)-P1_NEW(2))/(P4_NEW(3)-P1_NEW(3))

	PC(3)=0.D0

	R0=DSQRT(PC(1)**2+PC(2)**2+PC(3)**2)

	PC(:)=PC(:)*R/R0

	PCELL(1)=XX(1)*PC(1)+YY(1)*PC(2)
	PCELL(2)=XX(2)*PC(1)+YY(2)*PC(2)
	PCELL(3)=XX(3)*PC(1)+YY(3)*PC(2)

END SUBROUTINE CHANGE_CELL
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
SUBROUTINE FIND_LINE_NDIRECTION(PA,PB,PC)
		IMPLICIT NONE
		DOUBLE PRECISION,DIMENSION(1:3)::PA,PB,PC
		REAL*8 R,RC

		R=DSQRT(PA(1)**2+PA(2)**2+PA(3)**2)
		PC(1)=(PA(2)*PB(3)-PA(3)*PB(2))/R
		PC(2)=(PA(3)*PB(1)-PA(1)*PB(3))/R
		PC(3)=(PA(1)*PB(2)-PA(2)*PB(1))/R
		RC=DSQRT(PC(1)**2+PC(2)**2+PC(3)**2)
		IF(RC.GE.0.1D0**10)THEN
			PC=PC/RC
		ENDIF
ENDSUBROUTINE FIND_LINE_NDIRECTION
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
SUBROUTINE FIND_NEW_VECTOR(XX,YY,ZZ,P,PN)
	IMPLICIT NONE
	DOUBLE PRECISION,DIMENSION(1:3)::XX,YY,ZZ,P,PN

	PN(1)=XX(1)*P(1)+XX(2)*P(2)+XX(3)*P(3)
	PN(2)=YY(1)*P(1)+YY(2)*P(2)+YY(3)*P(3)
	PN(3)=ZZ(1)*P(1)+ZZ(2)*P(2)+ZZ(3)*P(3)
ENDSUBROUTINE FIND_NEW_VECTOR
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
SUBROUTINE FIND_LINE_DIRECTION(PA,PB,PC)	
		IMPLICIT NONE
		!PA,PB不重合 
		DOUBLE PRECISION,DIMENSION(1:3)::PA,PB,PC
		REAL*8 SITA,RA,RB
		REAL*8 T,A1,A2,A3,B1,B2,B3,FINE,R,TEST

		RA=DSQRT(PA(1)**2+PA(2)**2+PA(3)**2)
		RB=DSQRT(PB(1)**2+PB(2)**2+PB(3)**2)
		TEST=(PA(1)*PB(1)+PA(2)*PB(2)+PA(3)*PB(3))/RA/RB
		IF(TEST.GE.1.D0)THEN
			TEST=1.D0
		ENDIF
		IF(TEST.LE.-1.D0)THEN
			TEST=-1.D0
		ENDIF
		SITA=DACOS(TEST)

		A1=PA(1)/RA
		A2=PA(2)/RA
		A3=PA(3)/RA

		B1=PB(1)/RB
		B2=PB(2)/RB
		B3=PB(3)/RB

		T=DSIN(SITA)
		FINE=a2**2*b1**2 + a3**2*b1**2 - 2.D0*a1*a2*b1*b2 + a1**2*b2**2 + a3**2*b2**2 - 2.D0*a1*a3*b1*b3	&
		&		 - 2.D0*a2*a3*b2*b3 + a1**2*b3**2 + a2**2*b3**2

		PC(1)=(a2**2*b1*T + a3**2*b1*T - a1*a2*b2*T - a1*a3*b3*T)/FINE
		PC(2)=(-(a1*a2*b1*T) + a1**2*b2*T + a3**2*b2*T - a2*a3*b3*T)/FINE
		PC(3)=((-(a1*a3*b1) - a2*a3*b2 + a1**2*b3 + a2**2*b3)*T)/FINE

		R=DSQRT(PC(1)**2+PC(2)**2+PC(3)**2)
		PC(:)=PC(:)/R
ENDSUBROUTINE FIND_LINE_DIRECTION
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	SUBROUTINE INT_CHAR(N_OLD,CH)
		IMPLICIT NONE
		INTEGER(4) I,N,N0,P,P0,TEST,K,NN,N_OLD
		CHARACTER(20) CH
		CHARACTER(9)FILE

		N=N_OLD
		FILE='_mesh.dat'
		P=8
		N0=10**P

		K=0
                P0=0
		DO I=1,P
			TEST=INT(N/N0)
			IF(TEST==0)THEN
				N0=N0/10
				K=K+1
			ELSE
				P0=P-K
				EXIT
			ENDIF
		ENDDO

		NN=10**P0
		DO I=0,P0
			TEST=N/NN
			N=N-TEST*NN
			NN=NN/10
			CH(1+I:1+I)=CHAR(TEST+48)
		ENDDO
		CH(2+P0:10+P0)=FILE
	ENDSUBROUTINE INT_CHAR
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	SUBROUTINE CHANGE_NODECELLS(P,NODE,CELL,LINE_NODES,NODE_CELLS_OLD,NODE_LINES_OLD,CELL_LINES,NODE_CELLS,NODE_LINES,NODE_NODES)
		IMPLICIT NONE
		INTEGER P,I,J
		DOUBLE PRECISION,DIMENSION(1:20*P*P,1:3)::NODE
		DOUBLE PRECISION,DIMENSION(1:(12+30*(P-1)+(P-1)*(P-2)*10),1:3)::CELL
		INTEGER ,DIMENSION(1:(12+30*(P-1)+(P-1)*(P-2)*10),1:6)::CELL_NODES,CELL_LINES
		INTEGER ,DIMENSION(1:20*P*P,1:3)::NODE_CELLS,NODE_CELLS_OLD,NODE_LINES,NODE_LINES_OLD,NODE_NODES
		INTEGER ,DIMENSION(1:30*P*P,1:2)::LINE_NODES
		REAL*8 FINE
		INTEGER TEST

		DO I=1,20*P*P
			NODE_CELLS(I,1)=NODE_CELLS_OLD(I,1)
			CALL FIND_SECONDCELL(CELL(NODE_CELLS_OLD(I,1),:),NODE(I,:),CELL(NODE_CELLS_OLD(I,2),:),FINE)

			IF(FINE.GE.0.D0)THEN
				NODE_CELLS(I,2)=NODE_CELLS_OLD(I,2)
				NODE_CELLS(I,3)=NODE_CELLS_OLD(I,3)
			ELSE
				NODE_CELLS(I,2)=NODE_CELLS_OLD(I,3)
				NODE_CELLS(I,3)=NODE_CELLS_OLD(I,2)
			ENDIF

			TEST=(NODE_CELLS(I,1)-NODE_CELLS(I,2))
			IF(TEST.NE.0)TEST=1
			TEST=TEST*(NODE_CELLS(I,1)-NODE_CELLS(I,3))
			IF(TEST.NE.0)TEST=1
			TEST=TEST*(NODE_CELLS(I,2)-NODE_CELLS(I,3))
			IF(TEST==0)THEN
				WRITE(*,*)'CHANGE_NODECELLS ERROR'
				STOP
			ENDIF
		ENDDO

		DO I=1,20*P*P
			DO J=1,3
				CALL FIND_LINES(NODE_LINES_OLD(I,J),CELL_LINES(NODE_CELLS(I,1),:),CELL_LINES(NODE_CELLS(I,2),:),TEST)
				IF(TEST==0)THEN
					NODE_LINES(I,1)=NODE_LINES_OLD(I,J)
					EXIT
				ENDIF
			ENDDO
			DO J=1,3
				CALL FIND_LINES(NODE_LINES_OLD(I,J),CELL_LINES(NODE_CELLS(I,2),:),CELL_LINES(NODE_CELLS(I,3),:),TEST)
				IF(TEST==0)THEN
					NODE_LINES(I,2)=NODE_LINES_OLD(I,J)
					EXIT
				ENDIF
			ENDDO
			DO J=1,3
				CALL FIND_LINES(NODE_LINES_OLD(I,J),CELL_LINES(NODE_CELLS(I,3),:),CELL_LINES(NODE_CELLS(I,1),:),TEST)
				IF(TEST==0)THEN
					NODE_LINES(I,3)=NODE_LINES_OLD(I,J)
					EXIT
				ENDIF
			ENDDO
			TEST=(NODE_LINES(I,1)-NODE_LINES(I,2))
			IF(TEST.NE.0)TEST=1
			TEST=TEST*(NODE_LINES(I,1)-NODE_LINES(I,3))
			IF(TEST.NE.0)TEST=1
			TEST=TEST*(NODE_LINES(I,2)-NODE_LINES(I,3))
			IF(TEST==0)THEN
				WRITE(*,*)'CHANGE_NODECELLS ERROR'
				STOP
			ENDIF
		ENDDO

		DO I=1,20*P*P
			DO J=1,3
				IF(I-LINE_NODES(NODE_LINES(I,J),1)==0)THEN
					NODE_NODES(I,J)=LINE_NODES(NODE_LINES(I,J),2)
				ELSE
					NODE_NODES(I,J)=LINE_NODES(NODE_LINES(I,J),1)
				ENDIF
			ENDDO
		ENDDO
	ENDSUBROUTINE CHANGE_NODECELLS
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	SUBROUTINE FIND_SECONDCELL(C1,N0,C2,FINE)
		IMPLICIT NONE
		INTEGER I
		DOUBLE PRECISION,DIMENSION(1:3)::C1,N0,C2,NA,NB,NC
		REAL*8 FINE,R

		R=DSQRT(N0(1)**2+N0(2)**2+N0(3)**2)
		DO I=1,3
			NA(I)=(C1(I)-N0(I))/R
			NB(I)=(C2(I)-N0(I))/R
			NC(I)=N0(I)/R
		ENDDO
		FINE=(NA(2)*NB(3)-NA(3)*NB(2))*NC(1)+(NA(3)*NB(1)-NA(1)*NB(3))*NC(2)+(NA(1)*NB(2)-NA(2)*NB(1))*NC(3)
	ENDSUBROUTINE FIND_SECONDCELL
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	SUBROUTINE FIND_LINES(N0,CN1,CN2,TEST)
		IMPLICIT NONE
		INTEGER TEST,N0,I,T1,T2
		INTEGER,DIMENSION(1:6)::CN1,CN2

		T1=1
		DO I=1,6
			T1=T1*(N0-CN1(I))
			IF(T1.NE.0)THEN
				T1=1
			ELSE
				EXIT
			ENDIF
		ENDDO
		T2=1
		DO I=1,6
			T2=T2*(N0-CN2(I))
			IF(T2.NE.0)THEN
				T2=1
			ELSE
				EXIT
			ENDIF
		ENDDO
		TEST=T1+T2
	ENDSUBROUTINE FIND_LINES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	SUBROUTINE MESH_CENTER_POINT(P,NODE,LINE_NODES,NODE_ONLINE)
		IMPLICIT NONE
		INTEGER P,I
		DOUBLE PRECISION,DIMENSION(1:20*P*P,1:3)::NODE
		DOUBLE PRECISION,DIMENSION(1:30*P*P,1:3)::NODE_ONLINE

		INTEGER ,DIMENSION(1:30*P*P,1:2)::LINE_NODES

		DO I=1,30*P*P
			CALL FIND_CENTER_POINT(NODE(LINE_NODES(I,1),:),NODE(LINE_NODES(I,2),:),NODE_ONLINE(I,:))
		ENDDO
	ENDSUBROUTINE MESH_CENTER_POINT
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	SUBROUTINE FIND_CENTER_POINT(PA,PB,PC)
		IMPLICIT NONE
		INTEGER P,I
		DOUBLE PRECISION,DIMENSION(1:3):: PA,PB,PC
		REAL*8 SINTA
		REAL*8 KSAI,A,B,R
				
		SINTA=DACOS((PA(1)*PB(1)+PA(2)*PB(2)+PA(3)*PB(3)))
		KSAI=(PA(2)*PB(1)-PA(1)*PB(2))**2+(PA(2)*PB(3)-PA(3)*PB(2))**2+(PA(1)*PB(3)-PA(3)*PB(1))**2

		P=2
		I=2
		PC(1)=(DCOS(SINTA*DBLE(0.5))*(PA(1)*(PB(2)**2+PB(3)**2)-PB(1)*(PA(2)*PB(2)+PA(3)*PB(3)))+	&
			&	  DCOS(SINTA*DBLE(0.5))*(PB(1)*(PA(2)**2+PA(3)**2)-PA(1)*(PA(2)*PB(2)+PA(3)*PB(3))))/KSAI
		PC(2)=(DCOS(SINTA*DBLE(0.5))*(PA(2)*(PB(1)**2+PB(3)**2)-PB(2)*(PA(1)*PB(1)+PA(3)*PB(3)))+	&
			&	  DCOS(SINTA*DBLE(0.5))*(PB(2)*(PA(1)**2+PA(3)**2)-PA(2)*(PA(1)*PB(1)+PA(3)*PB(3))))/KSAI
		PC(3)=(DCOS(SINTA*DBLE(0.5))*(PA(3)*(PB(2)**2+PB(1)**2)-PB(3)*(PA(2)*PB(2)+PA(1)*PB(1)))+	&
			&	  DCOS(SINTA*DBLE(0.5))*(PB(3)*(PA(2)**2+PA(1)**2)-PA(3)*(PA(2)*PB(2)+PA(1)*PB(1))))/KSAI		
	ENDSUBROUTINE FIND_CENTER_POINT
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	SUBROUTINE MESH_LINES(P,NODE_CELLS,CELL_NODES,LINE_NODES,LINE_CELLS,CELL_LINES,NODE_LINES)
		IMPLICIT NONE
		INTEGER P,I,J,M,JJ,K,KK,L_N,LINE_N,FINE
		INTEGER ,DIMENSION(1:(12+30*(P-1)+(P-1)*(P-2)*10),1:6)::CELL_NODES
		INTEGER ,DIMENSION(1:20*P*P,1:3)::NODE_CELLS,NODE_LINES
		INTEGER ,DIMENSION(1:30*P*P,1:2)::LINE_NODES,LINE_CELLS
		INTEGER ,DIMENSION(1:(12+30*(P-1)+(P-1)*(P-2)*10),1:6)::CELL_LINES
		INTEGER ,DIMENSION(1:6*(12+30*(P-1)+(P-1)*(P-2)*10),1:2)::LINE_TEST

		DO I=1,(12+30*(P-1)+(P-1)*(P-2)*10)
			DO J=1,6
				IF(J==6)THEN
					JJ=1
				ELSE
					JJ=J+1
				ENDIF
				LINE_TEST((I-1)*6+J,1)=CELL_NODES(I,J)
				LINE_TEST((I-1)*6+J,2)=CELL_NODES(I,JJ)
			ENDDO
		ENDDO

		K=1
		DO I=1,(12+30*(P-1)+(P-1)*(P-2)*10)			
			DO J=1,6
				IF((LINE_TEST((I-1)*6+J,1)-LINE_TEST((I-1)*6+J,2)).NE.0)THEN
					L_N=0
					DO M=1,K
						CALL FIND_NEWLINE(LINE_TEST((I-1)*6+J,:),LINE_NODES(M,:),FINE)
						L_N=FINE
						IF(FINE.EQ.0)THEN
							LINE_N=M
							GOTO 2
						ENDIF
					ENDDO

2					IF(L_N.NE.0)THEN
						LINE_NODES(K,:)=LINE_TEST((I-1)*6+J,:)
						CELL_LINES(I,J)=K
						K=K+1
					ELSE
						CELL_LINES(I,J)=LINE_N
					ENDIF
				ENDIF						
			ENDDO	!J			
		ENDDO	!I

		DO I=1,20*P*P
			K=1
			DO J=1,30*P*P
				IF((LINE_NODES(J,1)-I)*(LINE_NODES(J,2)-I)==0)THEN
					NODE_LINES(I,K)=J
					K=K+1
				ENDIF
			ENDDO
			if(k.ne.4)then
				write(*,*)'Error'
				Stop
			endif
		ENDDO

		DO I=1,30*P*P
			JJ=1
			DO J=1,12+30*(P-1)+(P-1)*(P-2)*10
				DO K=1,6
					IF(I-CELL_LINES(J,K)==0)THEN
						LINE_CELLS(I,JJ)=J
						JJ=JJ+1
					ENDIF
				ENDDO
			!	IF(JJ==3) EXIT				
			ENDDO
			if(jj.ne.3)then
				write(*,*)'Error'
				Stop
			endif
		ENDDO
		ENDSUBROUTINE MESH_LINES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
		SUBROUTINE FIND_NEWLINE(L1,L2,FINE)
		IMPLICIT NONE 
		INTEGER FINE
		INTEGER,DIMENSION(1:2)::L1,L2

		FINE=ABS((L1(1)-L2(1))*(L1(1)-L2(2)))+ABS((L1(2)-L2(1))*(L1(2)-L2(2)))

		ENDSUBROUTINE FIND_NEWLINE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
		SUBROUTINE MESH(P,P_NODE_FINAL,CELL_NEW,NODE_CELLS_NEW,CELL_NODES_NEW)
		IMPLICIT NONE
		INTEGER P,I,J,K,L
		DOUBLE PRECISION,DIMENSION(1:3)::PNODE
		DOUBLE PRECISION,DIMENSION(1:12,1:3):: P_BASE
		DOUBLE PRECISION,DIMENSION(1:20,1:(P+1)*(P+2)/2,1:3):: P_CELL
		DOUBLE PRECISION,DIMENSION(1:20,1:P*P,1:3):: P_NODE
		DOUBLE PRECISION,DIMENSION(1:20*P*P,1:3):: P_NODE_FINAL		!FINE
		DOUBLE PRECISION,DIMENSION(1:10*(P+1)*(P+2),1:3):: P_CELL_FINAL  !FINE
		INTEGER,DIMENSION(1:20,1:P*P,1:3)::NODE_CELLS
		INTEGER,DIMENSION(1:20*P*P,1:3)::NODE_CELLS_FINAL			!FINE
		INTEGER,DIMENSION(1:20,1:(P+1)*(P+2)/2,1:6)::CELL_NODES
		INTEGER,DIMENSION(1:10*(P+1)*(P+2),1:6)::CELL_NODES_FINAL	!FINE
		INTEGER,DIMENSION(1:6)::S
		INTEGER,DIMENSION(1:9)::A
		INTEGER,DIMENSION(1:5)::B1_LINE,B2_LINE
		INTEGER,DIMENSION(1:5,1:3)::RIGHT_P,LEFT_P
		INTEGER,DIMENSION(6:15,1:3)::LINE3,PP
		INTEGER,DIMENSION(1:20,1:3)::CELL_BASE
		REAL*8 FINE
		REAL*8 RADIUS

		DOUBLE PRECISION,DIMENSION(1:20*P*P,1:3)::NODE_NEW
		DOUBLE PRECISION,DIMENSION(1:(12+30*(P-1)+(P-1)*(P-2)*10),1:3)::CELL_NEW
		INTEGER ,DIMENSION(1:(12+30*(P-1)+(P-1)*(P-2)*10),1:6)::CELL_NODES_NEW
		INTEGER ,DIMENSION(1:20*P*P,1:3)::NODE_CELLS_NEW
		INTEGER NEW_N,NEWCELL

		RADIUS=1.D0  !UNIT METER

		FINE=1.D0
		J=2
		K=3
		DO I=1,5
			CELL_BASE(I,1)=1
			CELL_BASE(I,2)=J
			CELL_BASE(I,3)=K
			J=J+1
			K=K+1
			IF(K.GT.6)K=2
		ENDDO
		J=8
		K=7
		DO I=1,5
			CELL_BASE(I+15,1)=7
			CELL_BASE(I+15,2)=J
			CELL_BASE(I+15,3)=K
			J=J+1
			K=K+1
			IF(I.EQ.1)CELL_BASE(I+15,3)=12
		ENDDO

		L=2
		J=3
		K=8
		DO I=6,14,2
			CELL_BASE(I,1)=L
			CELL_BASE(I,2)=J
			CELL_BASE(I,3)=K
			L=L+1
			J=J+1
			K=K+1
			IF(I.EQ.14)CELL_BASE(I,2)=2
		ENDDO

		L=3
		J=8
		K=9
		DO I=7,15,2
			CELL_BASE(I,1)=L
			CELL_BASE(I,2)=J
			CELL_BASE(I,3)=K
			L=L+1
			J=J+1
			K=K+1
			IF(I.EQ.15)THEN
				CELL_BASE(I,1)=2
				CELL_BASE(I,3)=8
			ENDIF
		ENDDO

		CALL BASE_POINT(P_BASE)
		
		CALL DIVIDE_TRIANGLE(P,P_BASE(1,:),P_BASE(2,:),P_BASE(3,:),P_CELL(1,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(1,:),P_BASE(3,:),P_BASE(4,:),P_CELL(2,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(1,:),P_BASE(4,:),P_BASE(5,:),P_CELL(3,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(1,:),P_BASE(5,:),P_BASE(6,:),P_CELL(4,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(1,:),P_BASE(6,:),P_BASE(2,:),P_CELL(5,:,:))

		CALL DIVIDE_TRIANGLE(P,P_BASE(2,:),P_BASE(12,:),P_BASE(8,:),P_CELL(6,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(8,:),P_BASE(3,:),P_BASE(2,:),P_CELL(7,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(3,:),P_BASE(8,:),P_BASE(9,:),P_CELL(8,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(9,:),P_BASE(4,:),P_BASE(3,:),P_CELL(9,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(4,:),P_BASE(9,:),P_BASE(10,:),P_CELL(10,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(10,:),P_BASE(5,:),P_BASE(4,:),P_CELL(11,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(5,:),P_BASE(10,:),P_BASE(11,:),P_CELL(12,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(11,:),P_BASE(6,:),P_BASE(5,:),P_CELL(13,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(6,:),P_BASE(11,:),P_BASE(12,:),P_CELL(14,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(12,:),P_BASE(2,:),P_BASE(6,:),P_CELL(15,:,:))

		CALL DIVIDE_TRIANGLE(P,P_BASE(7,:),P_BASE(8,:),P_BASE(12,:),P_CELL(16,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(7,:),P_BASE(9,:),P_BASE(8,:),P_CELL(17,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(7,:),P_BASE(10,:),P_BASE(9,:),P_CELL(18,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(7,:),P_BASE(11,:),P_BASE(10,:),P_CELL(19,:,:))
		CALL DIVIDE_TRIANGLE(P,P_BASE(7,:),P_BASE(12,:),P_BASE(11,:),P_CELL(20,:,:))

		DO L=1,20
			DO I=1,P
				K=(I+1)*(I+2)/2-I+1
				DO J=(I-1)**2+1,I**2,2
					CALL FIND_PNODE(P_CELL(L,K,:),P_CELL(L,K-1,:),P_CELL(L,K-(I+1),:),PNODE)
					P_NODE(L,J,:)=PNODE(:)
					NODE_CELLS(L,J,1)=K-1+(L-1)*(P+1)*(P+2)/2
					NODE_CELLS(L,J,2)=K+(L-1)*(P+1)*(P+2)/2
					NODE_CELLS(L,J,3)=K-(I+1)+(L-1)*(P+1)*(P+2)/2
					K=K+1
				ENDDO
			ENDDO
		ENDDO

		DO L=1,20
			DO I=2,P
				K=(I+1)*(I+2)/2-I+1
				DO J=(I-1)**2+2,I**2-1,2
					CALL FIND_PNODE(P_CELL(L,K,:),P_CELL(L,K-I,:),P_CELL(L,K-(I+1),:),PNODE)
					P_NODE(L,J,:)=PNODE(:)
					NODE_CELLS(L,J,1)=K+(L-1)*(P+1)*(P+2)/2
					NODE_CELLS(L,J,2)=K-I+(L-1)*(P+1)*(P+2)/2
					NODE_CELLS(L,J,3)=K-(I+1)+(L-1)*(P+1)*(P+2)/2
					K=K+1
				ENDDO
			ENDDO
		ENDDO

		DO I=1,20
			DO J=1,P*P
				P_NODE_FINAL((I-1)*P*P+J,:)=P_NODE(I,J,:)
				NODE_CELLS_FINAL((I-1)*P*P+J,:)=NODE_CELLS(I,J,:)
			ENDDO
		ENDDO

		DO I=1,20
			DO J=1,(P+1)*(P+2)/2
				P_CELL_FINAL((I-1)*(P+1)*(P+2)/2+J,:)=P_CELL(I,J,:)
			ENDDO
		ENDDO
!***************************************************************************
!FACE 1~5
		DO I=1,5
			B1_LINE(I)=7+2*(I-1)
			DO J=1,3
				IF(I.LE.4)THEN
					RIGHT_P(I,J)=7+2*(I-1)+(J-1)
				ELSE
					RIGHT_P(I,1)=15
					RIGHT_P(I,2)=6
					RIGHT_P(I,3)=7
				ENDIF
			ENDDO
		ENDDO

		DO L=1,5
			DO I=1,5
				CELL_NODES(L,1,I)=1+(I-1)*P*P
			ENDDO
			CELL_NODES(L,1,6)=CELL_NODES(L,1,1)
		ENDDO

		DO L=1,5
			K=(P+1)*(P+2)/2
			CELL_NODES(L,K,1)=L*P*P
			CELL_NODES(L,K,2)=(RIGHT_P(L,1)-1)*P*P+(P-1)**2+1
			CELL_NODES(L,K,3)=(RIGHT_P(L,2)-1)*P*P+1		
			CELL_NODES(L,K,4)=RIGHT_P(L,3)*P*P
			IF(L.LE.4)THEN
				CELL_NODES(L,K,5)=L*P*P+(P-1)**2+1
			ELSE
			!	CELL_NODES(L,K,5)=P*P+(P-1)**2+1
				CELL_NODES(L,K,5)=(P-1)**2+1
			ENDIF
			CELL_NODES(L,K,6)=CELL_NODES(L,K,1)
		ENDDO

		DO L=1,5
			K=(P+1)*(P+2)/2-P
			IF(L.EQ.1)THEN
				CELL_NODES(1,K,1)=(P-1)**2+1
				CELL_NODES(1,K,2)=5*P*P
				CELL_NODES(1,K,3)=14*P*P+1+(P-1)**2
				CELL_NODES(1,K,4)=5*P*P+1
				CELL_NODES(1,K,5)=7*P*P
				CELL_NODES(1,K,6)=CELL_NODES(1,K,1)
			ELSE
				CELL_NODES(L,K,:)=CELL_NODES(L-1,K+P,:)
			ENDIF
		ENDDO

	IF(P.GE.2)THEN
		DO L=1,5
!LINE 1
			DO I=1,P-1
				K=(I+1)*(I+2)/2-I
				IF(L.EQ.1)THEN
					CELL_NODES(1,K,1)=I*I+1
					CELL_NODES(1,K,2)=I*I+2
					CELL_NODES(1,K,3)=(I-1)**2+1
					CELL_NODES(1,K,4)=I*I+4*P*P
					CELL_NODES(1,K,5)=(I+1)**2-1+4*P*P
					CELL_NODES(1,K,6)=(I+1)**2+4*P*P
				ELSE
					CELL_NODES(L,K,:)=CELL_NODES(L-1,K+I,:)
				ENDIF
			ENDDO
!LINE 2
			DO I=1,P-1
				K=(I+1)*(I+2)/2				
				CELL_NODES(L,K,1)=(L-1)*P*P+I**2
				CELL_NODES(L,K,2)=(L-1)*P*P+(I+1)**2-1
				CELL_NODES(L,K,3)=(L-1)*P*P+(I+1)**2
				IF(L.LE.4)THEN
					CELL_NODES(L,K,4)=L*P*P+I**2+1
					CELL_NODES(L,K,5)=L*P*P+I**2+2
					CELL_NODES(L,K,6)=L*P*P+(I-1)**2+1
				ELSE
					CELL_NODES(L,K,4)=I*I+1
					CELL_NODES(L,K,5)=I*I+2
					CELL_NODES(L,K,6)=(I-1)**2+1
				ENDIF
			ENDDO
			
!LINE 3		
			J=0
			DO I=1,P-1
				K=(P+1)*(P+2)/2-P+I
				CELL_NODES(L,K,1)=(L-1)*P*P+(P-1)**2+3+J
				CELL_NODES(L,K,2)=(L-1)*P*P+(P-1)**2+2+J
				CELL_NODES(L,K,3)=(L-1)*P*P+(P-1)**2+1+J
				CELL_NODES(L,K,4)=(B1_LINE(L)-1)*P**2+P**2-J
				CELL_NODES(L,K,5)=(B1_LINE(L)-1)*P**2+P**2-1-J
				CELL_NODES(L,K,6)=(B1_LINE(L)-1)*P**2+P**2-2-J
				J=J+2
			ENDDO
		ENDDO
	ENDIF
!内点	
		IF(P.GE.3)THEN
			DO L=1,5
				DO I=2,P-1
					S(1)=(L-1)*P*P+(I-1)**2+3
					S(2)=(L-1)*P*P+(I-1)**2+2
					S(3)=(L-1)*P*P+(I-1)**2+1
					S(4)=(L-1)*P*P+I**2+2
					S(5)=(L-1)*P*P+I**2+3
					S(6)=(L-1)*P*P+I**2+4
					DO K=(I+1)*(I+2)/2-I+1,(I+1)*(I+2)/2-1
						CELL_NODES(L,K,:)=S(:)
						S(:)=S(:)+2
					ENDDO
				ENDDO				
			ENDDO
		ENDIF
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!FACE 16~20
		DO I=1,5
			B2_LINE(I)=6+2*(I-1)
			DO J=1,3
				IF(I.LE.4)THEN
					LEFT_P(I,J)=6+2*(I-1)+(J-1)
				ELSE
					LEFT_P(I,1)=14
					LEFT_P(I,2)=15
					LEFT_P(I,3)=6
				ENDIF
			ENDDO
		ENDDO

		DO L=16,20
			DO I=1,5
				CELL_NODES(L,1,7-I)=1+(15+I-1)*P*P
			ENDDO
			CELL_NODES(L,1,1)=CELL_NODES(L,1,6)           !???
		ENDDO
!LEFT
		DO L=16,20
			K=(P+1)*(P+2)/2-P
			CELL_NODES(L,K,6)=LEFT_P(L-15,1)*P**2
			CELL_NODES(L,K,5)=(LEFT_P(L-15,2)-1)*P*P+1
			CELL_NODES(L,K,4)=(LEFT_P(L-15,3)-1)*P*P+1+(P-1)**2
			IF(L.LE.19)THEN
				CELL_NODES(L,K,3)=(L+1)*P*P
			ELSE
				CELL_NODES(L,K,3)=16*P*P
			ENDIF
			CELL_NODES(L,K,2)=(L-1)*P*P+1+(P-1)**2
			CELL_NODES(L,K,1)=CELL_NODES(L,K,6)		!1~6	=> 6~1
		ENDDO
!RIGHT
		DO L=16,20
			K=(P+1)*(P+2)/2
			IF(L.EQ.16)THEN
				CELL_NODES(L,K,6)=L*P*P
				CELL_NODES(L,K,5)=5*P*P+(P-1)**2+1
				CELL_NODES(L,K,4)=14*P*P+1		
				CELL_NODES(L,K,3)=14*P*P
				CELL_NODES(L,K,2)=19*P*P+(P-1)**2+1
				CELL_NODES(L,K,1)=CELL_NODES(L,K,6)		!1~6	=> 6~1
			ELSE
				CELL_NODES(L,K,:)=CELL_NODES(L-1,K-P,:)
			ENDIF		
		ENDDO

	IF(P.GE.2)THEN
		DO L=16,20
!LINE 1
			DO I=1,P-1
				K=(I+1)*(I+2)/2-I
					CELL_NODES(L,K,1)=(L-1)*P*P+I*I+1
					CELL_NODES(L,K,2)=(L-1)*P*P+I*I+2
					CELL_NODES(L,K,3)=(L-1)*P*P+(I-1)**2+1
				IF(L.LE.19)THEN
					CELL_NODES(L,K,4)=L*P*P+I*I
					CELL_NODES(L,K,5)=L*P*P+(I+1)**2-1
					CELL_NODES(L,K,6)=L*P*P+(I+1)**2
				ELSE
					CELL_NODES(L,K,4)=15*P*P+I*I
					CELL_NODES(L,K,5)=15*P*P+(I+1)**2-1
					CELL_NODES(L,K,6)=15*P*P+(I+1)**2
				ENDIF
			ENDDO
!LINE 2
			DO I=1,P-1
				K=(I+1)*(I+2)/2				
				IF(L.EQ.16)THEN
					CELL_NODES(L,K,1)=(L-1)*P*P+I**2
					CELL_NODES(L,K,2)=(L-1)*P*P+(I+1)**2-1
					CELL_NODES(L,K,3)=(L-1)*P*P+(I+1)**2
					CELL_NODES(L,K,4)=19*P*P+I**2+1
					CELL_NODES(L,K,5)=19*P*P+I**2+2
					CELL_NODES(L,K,6)=19*P*P+(I-1)**2+1
				ELSE
					CELL_NODES(L,K,:)=CELL_NODES(L-1,K-I,:)
				ENDIF
			ENDDO
			
!LINE 3		
			J=0
			DO I=1,P-1
				K=(P+1)*(P+2)/2-P+I
				CELL_NODES(L,K,1)=(L-1)*P*P+(P-1)**2+3+J
				CELL_NODES(L,K,2)=(L-1)*P*P+(P-1)**2+2+J
				CELL_NODES(L,K,3)=(L-1)*P*P+(P-1)**2+1+J
				CELL_NODES(L,K,4)=(B2_LINE(L-15)-1)*P**2+P**2-J
				CELL_NODES(L,K,5)=(B2_LINE(L-15)-1)*P**2+P**2-1-J
				CELL_NODES(L,K,6)=(B2_LINE(L-15)-1)*P**2+P**2-2-J
				J=J+2
			ENDDO
		ENDDO
	ENDIF
!内点	
		IF(P.GE.3)THEN
			DO L=16,20
				DO I=2,P-1
					S(1)=(L-1)*P*P+(I-1)**2+3
					S(2)=(L-1)*P*P+(I-1)**2+2
					S(3)=(L-1)*P*P+(I-1)**2+1
					S(4)=(L-1)*P*P+I**2+2
					S(5)=(L-1)*P*P+I**2+3
					S(6)=(L-1)*P*P+I**2+4
					DO K=(I+1)*(I+2)/2-I+1,(I+1)*(I+2)/2-1
						CELL_NODES(L,K,:)=S(:)
						S(:)=S(:)+2
					ENDDO
				ENDDO				
			ENDDO
		ENDIF
!FACE 6~15	
		K=16
		DO L=6,14,2
			IF(L.EQ.6)THEN
				LINE3(L,1)=15
			ELSE
				LINE3(L,1)=L-1
			ENDIF
			LINE3(L,2)=L+1
			LINE3(L,3)=K
			K=K+1
		ENDDO	
		K=1
		DO L=7,15,2
			IF(L.EQ.15)THEN
				LINE3(L,1)=6
			ELSE
				LINE3(L,1)=L+1
			ENDIF
			LINE3(L,2)=L-1
			LINE3(L,3)=K
			K=K+1
		ENDDO	
!Vertices
		I=1
		J=15
		K=16
		DO L=6,14,2
			PP(L,1)=I
			IF(L.EQ.6)THEN
				PP(L,2)=20
			ELSE
				PP(L,2)=J
			ENDIF
			PP(L,3)=K
			I=I+1
			J=J+1
			K=K+1
		ENDDO
		I=16
		J=2
		K=1
		DO L=7,15,2
			PP(L,1)=I
			IF(L.EQ.15)THEN
				PP(L,2)=1
			ELSE
				PP(L,2)=J
			ENDIF
			PP(L,3)=K
			I=I+1
			J=J+1
			K=K+1
		ENDDO

		DO L=6,15
			CELL_NODES(L,1,:)=CELL_NODES(PP(L,1),(P+1)*(P+2)/2-P,:)
			CELL_NODES(L,(P+1)*(P+2)/2-P,:)=CELL_NODES(PP(L,2),(P+1)*(P+2)/2-P,:)
			CELL_NODES(L,(P+1)*(P+2)/2,:)=CELL_NODES(PP(L,3),(P+1)*(P+2)/2-P,:)
		ENDDO

		IF(P.GE.2)THEN
!LINE 1			
			DO L=6,15
				DO I=1,P-1
					K=(I+1)*(I+2)/2-I
					CELL_NODES(L,K,1)=(L-1)*P*P+I*I+1
					CELL_NODES(L,K,2)=(L-1)*P*P+I*I+2
					CELL_NODES(L,K,3)=(L-1)*P*P+(I-1)**2+1
					CELL_NODES(L,K,4)=(LINE3(L,1)-1)*P*P+(P-I)**2+1
					CELL_NODES(L,K,5)=(LINE3(L,1)-1)*P*P+(P-I)**2+2
					CELL_NODES(L,K,6)=(LINE3(L,1)-1)*P*P+(P-I-1)**2+1				
				ENDDO
!LINE 2
				DO I=1,P-1
					K=(I+1)*(I+2)/2				
					CELL_NODES(L,K,1)=(L-1)*P*P+I**2
					CELL_NODES(L,K,2)=(L-1)*P*P+(I+1)**2-1
					CELL_NODES(L,K,3)=(L-1)*P*P+(I+1)**2
					CELL_NODES(L,K,4)=(LINE3(L,2)-1)*P*P+(P-I)**2
					CELL_NODES(L,K,5)=(LINE3(L,2)-1)*P*P+(P-I+1)**2-1
					CELL_NODES(L,K,6)=(LINE3(L,2)-1)*P*P+(P-I+1)**2									
				ENDDO
!LINE 3			
				J=0
				DO I=1,P-1
					K=(P+1)*(P+2)/2-P+I
					CELL_NODES(L,K,1)=(L-1)*P*P+(P-1)**2+3+J
					CELL_NODES(L,K,2)=(L-1)*P*P+(P-1)**2+2+J
					CELL_NODES(L,K,3)=(L-1)*P*P+(P-1)**2+1+J
					CELL_NODES(L,K,4)=(LINE3(L,3)-1)*P**2+P**2-J
					CELL_NODES(L,K,5)=(LINE3(L,3)-1)*P**2+P**2-1-J
					CELL_NODES(L,K,6)=(LINE3(L,3)-1)*P**2+P**2-2-J
					J=J+2
				ENDDO
			ENDDO
		ENDIF

		IF(P.GE.3)THEN
			DO L=6,15
				DO I=2,P-1
					S(1)=(L-1)*P*P+(I-1)**2+3
					S(2)=(L-1)*P*P+(I-1)**2+2
					S(3)=(L-1)*P*P+(I-1)**2+1
					S(4)=(L-1)*P*P+I**2+2
					S(5)=(L-1)*P*P+I**2+3
					S(6)=(L-1)*P*P+I**2+4
					DO K=(I+1)*(I+2)/2-I+1,(I+1)*(I+2)/2-1
						CELL_NODES(L,K,:)=S(:)
						S(:)=S(:)+2
					ENDDO
				ENDDO				
			ENDDO
		ENDIF

		DO L=1,20
			DO I=1,(P+1)*(P+2)/2
				CELL_NODES_FINAL((L-1)*(P+1)*(P+2)/2+I,:)=CELL_NODES(L,I,:)
			ENDDO
		ENDDO

		DO I=1,20*P*P
			P_NODE_FINAL(I,:)=RADIUS*P_NODE_FINAL(I,:)
		ENDDO
		DO I=1,10*(P+1)*(P+2)
			P_CELL_FINAL(I,:)=RADIUS*P_CELL_FINAL(I,:)
		ENDDO

!------------------------------------------------CELL_NEW,CELL_NODES_NEW		
		K=(P+1)*(P+2)/2
		DO I=1,K
			CELL_NODES_NEW(I,:)=CELL_NODES_FINAL(I,:)
			CELL_NEW(I,:)=P_CELL_FINAL(I,:)
		ENDDO

		DO I=1+(P+1)*(P+2)/2,10*(P+1)*(P+2)	
			NEWCELL=0
			DO J=1,K
				CALL FIND_NEWCELL(CELL_NODES_FINAL(I,:),CELL_NODES_NEW(J,:),NEW_N)
				NEWCELL=NEW_N+NEWCELL
			ENDDO

			IF(NEWCELL.EQ.0)THEN ! (!=0) 表示重复
				K=K+1
				CELL_NODES_NEW(K,:)=CELL_NODES_FINAL(I,:)
				CELL_NEW(K,:)=P_CELL_FINAL(I,:)
			ENDIF
		ENDDO	
!------------------------------------------------NODE_CELLS_NEW
		DO I=1,20*P*P
			K=1
			DO J=1,12+30*(P-1)+(P-1)*(P-2)*10
				CALL FIND_NEWNODE(I,CELL_NODES_NEW(J,:),NEW_N)				
				IF(NEW_N==0)THEN					
					NODE_CELLS_NEW(I,K)=J
					K=K+1				
				ENDIF				
			ENDDO
			IF(K.NE.4)THEN
				WRITE(*,*)'	ERROR'
				WRITE(*,*)I,J,K
				STOP
			ENDIF
		ENDDO
!------------------------------------------------
!***************************************************************************		
		ENDSUBROUTINE MESH
!-----------------------------------------------------------------------------------
		SUBROUTINE FIND_NEWNODE(N0,NC,NEW_N)
		IMPLICIT NONE
		INTEGER NEW_N,N0,I
		INTEGER,DIMENSION(1:6)::NC

		NEW_N=1
		DO I=1,6
			NEW_N=NEW_N*ABS(N0-NC(I))
			IF(NEW_N.NE.0)THEN
				NEW_N=1
			ENDIF
		ENDDO

		ENDSUBROUTINE FIND_NEWNODE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		SUBROUTINE FIND_NEWCELL(CELL,CELL_NEW,NEW_N)
		IMPLICIT NONE
		INTEGER NEW_N,I,J,CT
		INTEGER,DIMENSION(1:6)::CELL,CELL_NEW

		NEW_N=0
		
		DO I=1,6
			CT=1			
			DO J=1,6
				CT=CT*ABS(CELL(I)-CELL_NEW(J))
				IF(CT.NE.0)THEN
					CT=1
				ENDIF			
			ENDDO
			IF(CT.NE.0) GOTO 1
		ENDDO
		NEW_N=1	!表示重复

1		ENDSUBROUTINE FIND_NEWCELL
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		SUBROUTINE DIVIDE_TRIANGLE(P,PA,PB,PC,PN)
		IMPLICIT NONE
		INTEGER P,I,J,K,KS,KE
		DOUBLE PRECISION,DIMENSION(1:3):: PA,PB,PC
		DOUBLE PRECISION,DIMENSION(1:(P+1)*(P+2)/2,1:3):: PN
		DOUBLE PRECISION,DIMENSION(1:(P+1),1:3):: PN_LINE
		REAL*8,ALLOCATABLE::PI_LINE(:,:)

		CALL DIVIDE_LINE(P,PA,PB,PN_LINE)	
		DO I=1,P+1
			K=(I+1)*I/2-(I-1)
			PN(K,:)=PN_LINE(I,:)			
		ENDDO

		CALL DIVIDE_LINE(P,PA,PC,PN_LINE)
		DO I=1,P+1
			K=(I+1)*I/2
			PN(K,:)=PN_LINE(I,:)			
		ENDDO

		DO I=1,P
			KS=(I+1)*(I+2)/2-I
			KE=(I+1)*(I+2)/2
			ALLOCATE(PI_LINE(1:I+1,1:3))
			CALL DIVIDE_LINE(I,PN(KS,:),PN(KE,:),PI_LINE)
			DO J=1,I+1
				PN((I+1)*(I+2)/2-I+J-1,:)=PI_LINE(J,:)
			ENDDO
			DEALLOCATE(PI_LINE)
		ENDDO
		
		ENDSUBROUTINE DIVIDE_TRIANGLE
!-----------------------------------------------------------------------------------
		SUBROUTINE DIVIDE_LINE(P,PA,PB,PN_LINE)
		IMPLICIT NONE
		INTEGER P,I
		DOUBLE PRECISION,DIMENSION(1:3):: PA,PB
		DOUBLE PRECISION,DIMENSION(1:(P+1),1:3):: PN_LINE
		REAL*8 SINTA
		REAL*8 KSAI,A,B

		SINTA=DACOS(PA(1)*PB(1)+PA(2)*PB(2)+PA(3)*PB(3))
		KSAI=(PA(2)*PB(1)-PA(1)*PB(2))**2+(PA(2)*PB(3)-PA(3)*PB(2))**2+(PA(1)*PB(3)-PA(3)*PB(1))**2
		PN_LINE(1,:)=PA
		PN_LINE(P+1,:)=PB
		IF(P.GE.2)THEN
			DO I=2,P
				PN_LINE(I,1)=(DCOS(SINTA*DBLE(1.D0*(I-1)/P))*(PA(1)*(PB(2)**2+PB(3)**2)-PB(1)*(PA(2)*PB(2)+PA(3)*PB(3)))+	&
					&	  DCOS(SINTA*DBLE(1.D0*(P-I+1)/P))*(PB(1)*(PA(2)**2+PA(3)**2)-PA(1)*(PA(2)*PB(2)+PA(3)*PB(3))))/KSAI
				PN_LINE(I,2)=(DCOS(SINTA*DBLE(1.D0*(I-1)/P))*(PA(2)*(PB(1)**2+PB(3)**2)-PB(2)*(PA(1)*PB(1)+PA(3)*PB(3)))+	&
					&	  DCOS(SINTA*DBLE(1.D0*(P-I+1)/P))*(PB(2)*(PA(1)**2+PA(3)**2)-PA(2)*(PA(1)*PB(1)+PA(3)*PB(3))))/KSAI
				PN_LINE(I,3)=(DCOS(SINTA*DBLE(1.D0*(I-1)/P))*(PA(3)*(PB(2)**2+PB(1)**2)-PB(3)*(PA(2)*PB(2)+PA(1)*PB(1)))+	&
					&	  DCOS(SINTA*DBLE(1.D0*(P-I+1)/P))*(PB(3)*(PA(2)**2+PA(1)**2)-PA(3)*(PA(2)*PB(2)+PA(1)*PB(1))))/KSAI
			ENDDO
		ENDIF
		
		ENDSUBROUTINE DIVIDE_LINE
!-----------------------------------------------------------------------------------
		SUBROUTINE BASE_POINT(PBASE)
		IMPLICIT NONE
		INTEGER I
		DOUBLE PRECISION,DIMENSION(1:12,1:3):: PBASE
		DOUBLE PRECISION,DIMENSION(1:12):: X,Y,Z
		REAL*8 A,B,C,D,E

		A=0.2D0*DSQRT(5.D0)
		B=(1.D0-A)*0.5D0
		C=(1.D0+A)*0.5D0
		D=DSQRT(B)
		E=DSQRT(C)

		X(1)=0.D0
		Y(1)=0.D0
		Z(1)=1.D0

		X(2)=2.D0*A
		Y(2)=0.D0
		Z(2)=A

		X(3)=B
		Y(3)=E
		Z(3)=A

		X(4)=-C
		Y(4)=D
		Z(4)=A

		X(5)=-C
		Y(5)=-D
		Z(5)=A

		X(6)=B
		Y(6)=-E
		Z(6)=A

		X(7)=0.D0
		Y(7)=0.D0
		Z(7)=-1.D0

		X(8)=C
		Y(8)=D
		Z(8)=-A

		X(9)=-B
		Y(9)=E
		Z(9)=-A

		X(10)=-2.D0*A
		Y(10)=0.D0
		Z(10)=-A

		X(11)=-B
		Y(11)=-E
		Z(11)=-A

		X(12)=C
		Y(12)=-D
		Z(12)=-A

		DO I=1,12
			PBASE(I,1)=X(I)
			PBASE(I,2)=Y(I)
			PBASE(I,3)=Z(I)
		ENDDO

		ENDSUBROUTINE BASE_POINT
!--------------------------------------------------------------
		SUBROUTINE FIND_PNODE(P1,P2,P3,PNODE)
		IMPLICIT NONE
		DOUBLE PRECISION,DIMENSION(1:3):: P1,P2,P3,PNODE,P_TEXT1,P_TEXT2
		DOUBLE PRECISION,DIMENSION(1:4,1:3):: P_LINE
		DOUBLE PRECISION,DIMENSION(1:3,1:3):: PT_LINE
		
		CALL DIVIDE_LINE(3,P1,P2,P_LINE)
		P_TEXT1(:)=P_LINE(3,:)

		CALL DIVIDE_LINE(3,P1,P3,P_LINE)
		P_TEXT2(:)=P_LINE(3,:)

		CALL DIVIDE_LINE(2,P_TEXT1,P_TEXT2,PT_LINE)

		PNODE(:)=PT_LINE(2,:)
		ENDSUBROUTINE FIND_PNODE
!--------------------------------------------------------------	
