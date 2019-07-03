#include "../hopest_f.h"

MODULE MODH_Mesh
!===================================================================================================================================
! Contains subroutines to build (curviilinear) meshes and provide metrics, etc.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitMesh
  MODULE PROCEDURE InitMesh
END INTERFACE

! INTERFACE SetCurvedInfo
!   MODULE PROCEDURE SetCurvedInfo
! END INTERFACE

! INTERFACE BuildHOMesh
!   MODULE PROCEDURE BuildHOMesh
! END INTERFACE

! INTERFACE DeformMesh
!   MODULE PROCEDURE DeformMesh
! END INTERFACE

INTERFACE FinalizeMesh
  MODULE PROCEDURE FinalizeMesh
END INTERFACE

PUBLIC::InitMesh
! PUBLIC::SetCurvedInfo
! PUBLIC::BuildHOMesh
! PUBLIC::DeformMesh
PUBLIC::FinalizeMesh
!===================================================================================================================================

CONTAINS


! Actually we don't need it, ProjectName and MeshFile have Already been defined 
SUBROUTINE InitMesh()
!===================================================================================================================================
! Read Parameter from inputfile 
!===================================================================================================================================
! MODULES
USE MOD_Globals
! USE MODH_Output_Vars, ONLY: ProjectName
! USE MODH_Mesh_Vars,   ONLY: Deform
! USE MODH_Mesh_Vars,   ONLY: MeshFile,Deform
! USE MODH_Mesh_Vars,   ONLY: doSplineInterpolation
! USE MODH_ReadInTools, ONLY: GETINT,CNTSTR,GETLOGICAL
USE MOD_ReadInTools, ONLY: GETSTR
USE MOD_Mesh_Vars,   ONLY: MeshFile
! USE MOD_Output_vars,       ONLY: ProjectName

IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT AMR MESH ...'
! prepare pointer structure (get nTrees, etc.)

MeshFile = GETSTR('MeshFile')

! Deform = GETINT('Deform','0')


! doSplineInterpolation = GETLOGICAL('doSplineInterpolation','.FALSE.')

SWRITE(UNIT_stdOut,'(A)')' INIT AMR MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitMesh


! SUBROUTINE SetCurvedInfo()
! !===================================================================================================================================
! ! Set and allocate information related to high order data
! !===================================================================================================================================
! ! MODULES
! USE MODH_Globals
! USE MOD_Mesh_Vars,ONLY: NGeo
! USE MODH_Mesh_Vars,ONLY: Xi_NGeo,wBary_NGeo,HexMap,HexMapInv
! USE MODH_Mesh_Vars,ONLY: NGeo_out,XiCL_NGeo_out,wBaryCL_Ngeo_out,HexMap_out
! USE MODH_Mesh_Vars,ONLY: Vdm_01,Vdm_10,Vdm_CL_EQ_out
! USE MODH_Mesh_Vars,ONLY: nCurvedNodes 
! USE MODH_Basis,    ONLY: BarycentricWeights,ChebyGaussLobNodesAndWeights,InitializeVandermonde
! USE MODH_ReadInTools, ONLY: GETINT

! !-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT NONE
! ! INPUT VARIABLES
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT/OUTPUT VARIABLES
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! OUTPUT VARIABLES
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! INTEGER :: i,j,k,l
! REAL,ALLOCATABLE    :: xi_Ngeo_out(:)
! CHARACTER(LEN=5) :: tmpstr
! !===================================================================================================================================
! ALLOCATE(Xi_Ngeo(0:NGeo))
! ALLOCATE(wBary_Ngeo(0:Ngeo))
! DO i=0,NGeo
!   Xi_Ngeo(i)=-1+REAL(i)*2./REAL(NGeo)
! END DO
! CALL BarycentricWeights(Ngeo,xi_Ngeo,wBary_Ngeo)

! !dirty readin
! ! WRITE(tmpstr,'(I5)')Ngeo
! ! Ngeo_out = GETINT('Ngeo_out',tmpstr)
! Ngeo_out = Ngeo

! Ngeo_out = MIN(Ngeo,Ngeo_out) ! should be at maximum Ngeo

! ALLOCATE(Xi_Ngeo_out(0:Ngeo_out),wBaryCL_Ngeo_out(0:Ngeo_out))
! ALLOCATE(XiCL_Ngeo_out(0:Ngeo_out))
! CALL ChebyGaussLobNodesAndWeights(Ngeo_out,XiCL_Ngeo_out)
! CALL BarycentricWeights(Ngeo_out,xiCL_Ngeo_out,wBaryCL_Ngeo_out)

! !for output from CL to Equidistant points
! DO i=0,NGeo_out
!   Xi_Ngeo_out(i)=-1+REAL(i)*2./REAL(NGeo_out)
! END DO

! !only used in output!!!
! ALLOCATE(Vdm_CL_EQ_out(0:Ngeo_out,0:Ngeo_out))
! CALL InitializeVandermonde(Ngeo_out,Ngeo_out,wBaryCL_Ngeo_out,xiCL_Ngeo_out,xi_Ngeo_Out,Vdm_CL_EQ_out)


! ALLOCATE(Vdm_10(0:Ngeo_out,0:NGeo_out),Vdm_01(0:Ngeo_out,0:NGeo_out))
! ! change from interval [-1,1] -> [-1,0] Vdm_10 and interval [-1,1]-> [0,1] Vdm__01
! CALL InitializeVandermonde(Ngeo_out,Ngeo_out,wBaryCL_Ngeo_out,xiCL_Ngeo_out,-1+0.5*(xiCL_Ngeo_Out+1),Vdm_10)
! CALL InitializeVandermonde(Ngeo_out,Ngeo_out,wBaryCL_Ngeo_out,xiCL_Ngeo_out,   0.5*(xiCL_Ngeo_Out+1),Vdm_01)

! ! mapping form one-dimensional list [1 ; (Ngeo+1)^3] to tensor-product 0 <= i,j,k <= Ngeo and back
! ALLOCATE(HexMap(0:Ngeo,0:Ngeo,0:Ngeo),HexMapInv(3,(Ngeo+1)**3))
! l=0
! DO k=0,Ngeo ; DO j=0,Ngeo ; DO i=0,Ngeo
!   l=l+1
!   HexMap(i,j,k)=l
!   HexMapInv(:,l)=(/i,j,k/)
! END DO ; END DO ; END DO
! ALLOCATE(HexMap_out(0:Ngeo_out,0:Ngeo_out,0:Ngeo_out))
! l=0
! DO k=0,Ngeo_out ; DO j=0,Ngeo_out ; DO i=0,Ngeo_out
!   l=l+1
!   HexMap_out(i,j,k)=l
! END DO ; END DO ; END DO

! IF(NGeo.GT.1)THEN
!   nCurvedNodes=(NGeo+1)**3
! ELSE
!   nCurvedNodes=0
! END IF

! END SUBROUTINE SetCurvedInfo


! SUBROUTINE BuildHOMesh()
! !===================================================================================================================================
! ! uses XGeo High order data from trees and interpolates it to the quadrants 
! !===================================================================================================================================
! ! MODULES
! USE MODH_Globals
! USE MODH_Mesh_Vars,   ONLY: Ngeo,nTrees,nElems,Xgeo,XgeoElem
! USE MODH_Mesh_Vars,   ONLY: nTrees,Trees
! USE MODH_Mesh_Vars,   ONLY: wBary_Ngeo,xi_Ngeo
! USE MODH_Mesh_Vars,   ONLY: wBaryCL_Ngeo_out
! USE MODH_Mesh_Vars,   ONLY: doSplineInterpolation
! USE MODH_Spline1D,    ONLY: GetSplineVandermonde
! USE MODH_P4EST_Vars,  ONLY: QuadToTree,QuadCoords,QuadLevel,sIntSize
! USE MODH_Basis,       ONLY: LagrangeInterpolationPolys,InitializeVandermonde
! USE MODH_ChangeBasis, ONLY: ChangeBasis3D, ChangeBasis3D_XYZ 
! USE MODH_ChangeBasis, ONLY: ChangeBasis2D_XY
! USE MODH_Mesh_Vars,   ONLY: Elems,BoundaryType
! USE MODH_Mesh_Vars,   ONLY: Ngeo_out,xiCL_Ngeo_out
! USE MODH_Mesh_Vars,   ONLY: blending_glob
! USE MODH_P4EST_Vars,  ONLY: P2H_FaceMap,P2H_VertexMap
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT VARIABLES
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! OUTPUT VARIABLES
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! INTEGER                           :: i,iElem,iTree 
! REAL                              :: xi0(3)
! REAL                              :: dxi,length
! REAL,DIMENSION(0:Ngeo_out,0:Ngeo) :: Vdm_Ngeo_NGeo_out
! REAL,DIMENSION(0:Ngeo_out,0:Ngeo_out) :: Vdm_xi,Vdm_eta,Vdm_zeta
! REAL                              :: XGeo_tree(3,0:NGeo_out,0:NGeo_out,0:NGeo_out)
! !-----------------------------------------------------------------------------------------------------------------------------------
! LOGICAL                           :: BCSide(0:5)
! INTEGER                           :: j,k,plus,BCIndex
! INTEGER                           :: PlocSide,PoppSide,nBCs,iNode,NodeSum 
! INTEGER                           :: MortarType(0:5)
! REAL,DIMENSION(0:Ngeo_out,0:Ngeo) :: Vdm_xi2,Vdm_eta2,Vdm_zeta2
! REAL                              :: NodeMark(0:1,0:1,0:1),lengthIJK(3),xi0IJK(3),xi(3)
! REAL                              :: xGeoElem_corr(3,0:NGeo_out,0:NGeo_out,0:NGeo_out)
! REAL                              :: blending(0:NGeo_out,0:NGeo_out,0:NGeo_out)
! !===================================================================================================================================
! ALLOCATE(XgeoElem(3,0:Ngeo_out,0:Ngeo_out,0:Ngeo_out,nElems))
! ALLOCATE(blending_glob(1,0:Ngeo_out,0:Ngeo_out,0:Ngeo_out,nElems))
! blending_glob=0.

! DO iElem=1,nElems
!   iTree=QuadToTree(iElem)+1
!   !interpolate tree HO mapping from NGeo to Ngeo_out
!   CALL InitializeVandermonde(Ngeo,Ngeo_out,wBary_Ngeo,xi_Ngeo,xiCL_Ngeo_out,Vdm_Ngeo_Ngeo_out)
!   CALL ChangeBasis3D(3,Ngeo,Ngeo_out,Vdm_Ngeo_Ngeo_out,XGeo(:,:,:,:,iTree),Xgeo_tree)

!   ! transform p4est first corner coordinates (integer from 0... intsize) to [-1,1] reference element
!   xi0(:)=-1.+2.*REAL(QuadCoords(:,iElem))*sIntSize
!   ! length of each quadrant in integers
!   length=2./REAL(2**QuadLevel(iElem))
!   ! Build Vandermonde matrices for each parameter range in xi, eta,zeta
!   IF(doSplineInterpolation)THEN
!     CALL getSplineVandermonde(Ngeo_out+1,Ngeo_out+1,Vdm_xi(:,:)  ,xi0(1)+0.5*(xiCL_Ngeo_out(:)+1.)*Length)
!     CALL getSplineVandermonde(Ngeo_out+1,Ngeo_out+1,Vdm_eta(:,:) ,xi0(2)+0.5*(xiCL_Ngeo_out(:)+1.)*Length)
!     CALL getSplineVandermonde(Ngeo_out+1,Ngeo_out+1,Vdm_zeta(:,:),xi0(3)+0.5*(xiCL_Ngeo_out(:)+1.)*Length)
!   ELSE !polynomials
!     DO i=0,Ngeo_out
!       dxi=0.5*(xiCL_Ngeo_out(i)+1.)*Length
!       CALL LagrangeInterpolationPolys(xi0(1) + dxi,Ngeo_out,xiCL_Ngeo_out,wBaryCL_Ngeo_out,Vdm_xi(i,:)) 
!       CALL LagrangeInterpolationPolys(xi0(2) + dxi,Ngeo_out,xiCL_Ngeo_out,wBaryCL_Ngeo_out,Vdm_eta(i,:)) 
!       CALL LagrangeInterpolationPolys(xi0(3) + dxi,Ngeo_out,xiCL_Ngeo_out,wBaryCL_Ngeo_out,Vdm_zeta(i,:)) 
!     END DO
!   END IF
!   !interpolate tree HO mapping to quadrant HO mapping (If Ngeo_out < Ngeo: Interpolation error!)
!   CALL ChangeBasis3D_XYZ(3,Ngeo_out,Ngeo_out,Vdm_xi,Vdm_eta,Vdm_zeta,XGeo_tree,XgeoElem(:,:,:,:,iElem))
! END DO !iElem=1,nElems

! IF(Ngeo_out.GE.Ngeo) THEN
!   RETURN
! ELSE
!  WRITE(*,*)'!!!!!!!!! WARNING: Correction of the BC surfaces with quad-local mapping!'
! END IF
! ! For Ngeo>1 and Ngeo_out < Ngeo, we need to correct the mortar faces!!

! ! correction of the BC faces (evaluation of the Ngeo Tree Mapping on Quads and watertight)

! !reset the node pointer tmp variable
! DO iTree=1,nTrees
!   DO iNode=1,8
!     Trees(iTree)%ep%Node(iNode)%np%tmp=0
!   END DO! iNode=1,8
! END DO! iTree=1,nTrees
! ! mark tree nodes for blending
! DO iTree=1,nTrees
!   DO PlocSide=0,5
!     BCIndex=Trees(iTree)%ep%Side(P2H_FaceMap(PlocSide))%sp%BCIndex
!     IF(BCIndex.GT.0) THEN
!       DO iNode=1,4
!         Trees(iTree)%ep%Side(P2H_FaceMap(PlocSide))%sp%Node(iNode)%np%tmp=1
!       END DO! iNode=1,4
!     END IF
!   END DO! PlocSide=0,5
! END DO! iTree=1,nTrees


! DO iElem=1,nElems
!   IF(QuadLevel(iElem).EQ.0) CYCLE !we dont correct level 0
!   IF(QuadLevel(iElem).GT.1) STOP 'correction only for level 1!!!'
!   iTree=QuadToTree(iElem)+1
!   nBCs=0
!   ! check if we need to correct the tree
!   NodeSum=0
!   DO iNode=1,8
!     NodeSum=NodeSum+Trees(iTree)%ep%Node(iNode)%np%tmp
!   END DO! iNode=1,8
!   IF(NodeSum.EQ.0) CYCLE
!   ! check if the tree has a BC side
!   BCSide=.FALSE.
!   DO PlocSide=0,5
!     BCIndex=Trees(iTree)%ep%Side(P2H_FaceMap(PlocSide))%sp%BCIndex
!     IF (BCIndex.NE.0) THEN 
!       IF (BoundaryType(1,BCIndex).NE.1) THEN ! we dont want to correct periodic BC
!         BCSide(PlocSide)=.TRUE.
!         nBCs=nBCs+1
!       END IF
!     END IF
!     MortarType(PlocSide)=Elems(iElem)%ep%Side(P2H_FaceMap(PLocSide))%sp%MortarType
!   END DO

!   ! interpolate the HO volume mapping from tree level+Ngeo to quad level+Ngeo_out
!   xi0(:)=-1.+2.*REAL(QuadCoords(:,iElem))*sIntSize
!   ! length of each quadrant in integers
!   length=2./REAL(2**QuadLevel(iElem))
!   ! Build Vandermonde matrices for each parameter range in xi, eta,zeta
!   IF(doSplineInterpolation)THEN
!     CALL getSplineVandermonde(Ngeo+1,Ngeo_out+1,Vdm_xi2(:,:)  ,xi0(1)+0.5*(xiCL_Ngeo_out(:)+1.)*Length)
!     CALL getSplineVandermonde(Ngeo+1,Ngeo_out+1,Vdm_eta2(:,:) ,xi0(2)+0.5*(xiCL_Ngeo_out(:)+1.)*Length)
!     CALL getSplineVandermonde(Ngeo+1,Ngeo_out+1,Vdm_zeta2(:,:),xi0(3)+0.5*(xiCL_Ngeo_out(:)+1.)*Length)
!   ELSE !polynomials
!     DO i=0,Ngeo_out
!       dxi=0.5*(xiCL_Ngeo_out(i)+1.)*Length
!       CALL LagrangeInterpolationPolys(xi0(1) + dxi,Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_xi2(i,:)) 
!       CALL LagrangeInterpolationPolys(xi0(2) + dxi,Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_eta2(i,:)) 
!       CALL LagrangeInterpolationPolys(xi0(3) + dxi,Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_zeta2(i,:)) 
!     END DO
!   END IF
!   !interpolate tree HO mapping on quadrant at Ngeo_out
!   CALL ChangeBasis3D_XYZ(3,Ngeo,Ngeo_out,Vdm_xi2,Vdm_eta2,Vdm_zeta2,XGeo(:,:,:,:,iTree),xGeoElem_corr)


!   ! generate the blending function based on the previously marked nodes and the mortar status of the faces
!   nodeMark=0.
!   lengthIJK=length
!   xi0IJK=xi0
!   iNode=0
!   DO k=0,1
!     DO j=0,1
!       DO i=0,1
!         nodeMark(i,j,k)= REAL(Trees(iTree)%ep%Node(P2H_VertexMap(iNode))%np%tmp)
!         iNode=iNode+1
!       END DO! i=0,1
!     END DO! j=0,1
!   END DO! k=0,1

!   DO PLocSide=0,5
!     plus=MOD(PlocSide,2) ! plus=0: minus side, plus=1: plus Side 
!     PoppSide=PlocSide +1 - 2*plus
!     IF((MortarType(PlocSide).LT.0) .AND. (.NOT.BCSide(PoppSide)))THEN
!       SELECT CASE(PlocSide)
!         CASE(0,1) !ximinus,xiplus
!           nodeMark(plus,:,:)=0.
!           lengthIJK(1)=2.
!           xi0IJK(1)   =-1.
!         CASE(2,3) !etaminus,etaplus
!           nodeMark(:,plus,:)=0.
!           lengthIJK(2)=2.
!           xi0IJK(2)   =-1.
!         CASE(4,5) !zetaminus,zetaplus
!           nodeMark(:,:,plus)=0.
!           lengthIJK(3)=2.
!           xi0IJK(3)   =-1.
!       END SELECT !PlocSide
!     END IF!(BCSide(PlocSide))
!   END DO! PLocSide=0,5

!   ! map the node values to the volume cell
!   ! xi is in [0,1] within the tree (BC dir) and in[0,1] within the element (smallMortarDir...)
!   DO k=0,NGeo_out
!     xi(3)=0.5*(1.+ xi0IJK(3) +lengthIJK(3)*0.5*(1.+xiCL_NGeo_out(k)))
!     DO j=0,NGeo_out
!       xi(2)=0.5*(1.+ xi0IJK(2) +lengthIJK(2)*0.5*(1.+xiCL_NGeo_out(j)))
!       DO i=0,NGeo_out
!         xi(1)=0.5*(1.+ xi0IJK(1) +lengthIJK(1)*0.5*(1.+xiCL_NGeo_out(i)))
!         blending(i,j,k)=nodeMark(0,0,0)*(1.-xi(1))*(1.-xi(2))*(1.-xi(3))  &
!                        +nodeMark(1,0,0)*xi(1)     *(1.-xi(2))*(1.-xi(3))  &
!                        +nodeMark(0,1,0)*(1.-xi(1))*xi(2)     *(1.-xi(3))  &
!                        +nodeMark(1,1,0)*xi(1)     *xi(2)     *(1.-xi(3))  &
!                        +nodeMark(0,0,1)*(1.-xi(1))*(1.-xi(2))*xi(3)       &
!                        +nodeMark(1,0,1)*xi(1)     *(1.-xi(2))*xi(3)       &
!                        +nodeMark(0,1,1)*(1.-xi(1))*xi(2)     *xi(3)       &
!                        +nodeMark(1,1,1)*xi(1)     *xi(2)     *xi(3)        
!       END DO! i=0,NGeo_out
!     END DO! j=0,NGeo_out
!   END DO! k=0,NGeo_out
!   blending_glob(1,:,:,:,iElem)=blending
!   ! blend the HO mapping with the old LO mapping
!   DO i=1,3 
!     xGeoElem(i,:,:,:,iElem)=xGeoElem(i,:,:,:,iElem)*(1.-blending(:,:,:))&
!                            +xGeoElem_corr(i,:,:,:)*blending(:,:,:)
!   END DO! i=1,3 
! END DO !iElem=1,nElems

! END SUBROUTINE BuildHOMesh



! SUBROUTINE TransFace(dim1,N_in,xi_in,Face)
! !===================================================================================================================================
! ! Transfinite mapping edge faces -> face, replace inner points only 
! !===================================================================================================================================
! ! MODULES
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT VARIABLES
! INTEGER,INTENT(IN)         :: dim1           ! size of leading dimension
! INTEGER,INTENT(IN)         :: N_in           !polynomial degree
! REAL,INTENT(IN)            :: xi_in(0:N_in)  !node positions in parameter space [-1,1]
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! OUTPUT VARIABLES
! REAL,INTENT(INOUT)         :: Face(1:dim1,0:N_in,0:N_in) !Face data, inner points will be overwritten
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! INTEGER                        :: i,j
! REAL                           :: xi0(0:N_in),xiN(0:N_in)
! !-----------------------------------------------------------------------------------------------------------------------------------
! IF(N_in.EQ.1) RETURN
! xiN=0.5*(xi_in+1.)  ![-1,1] ->  0...1
! xi0=1.-xiN          ! 1 ... 0
! DO j=1,N_in-1
!   DO i=1,N_in-1
!      Face(:,i,j)=   Face(:,i,0)*xi0(j) + Face(:,i,N_in)*xiN(j)  & !edges
!                    +Face(:,0,j)*xi0(i) + Face(:,N_in,j)*xiN(i)  & 
!                   -(  Face(:,   0,   0)*xi0(i)*xi0(j)          & !-corners
!                     + Face(:,N_in,   0)*xiN(i)*xi0(j)          &
!                     + Face(:,   0,N_in)*xi0(i)*xiN(j)          &
!                     + Face(:,N_in,N_in)*xiN(i)*xiN(j))
!   END DO !i
! END DO !j
! END SUBROUTINE TransFace


! SUBROUTINE TransVol(dim1,N_in,xi_in,Vol_in,Vol)
! !===================================================================================================================================
! ! Transfinite mapping 6 faces -> volume, replace face points and move inner points by a transfinite interpolation of face deformation
! !===================================================================================================================================
! ! MODULES
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT VARIABLES
! INTEGER,INTENT(IN)         :: dim1           ! size of leading dimension
! INTEGER,INTENT(IN)         :: N_in           !polynomial degree
! REAL,INTENT(IN)            :: xi_in(0:N_in)  !node positions in parameter space [-1,1]
! REAL,INTENT(IN)            :: Vol_in(1:dim1,0:N_in,0:N_in,0:N_in) ! volume data, but only face data is used
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! OUTPUT VARIABLES
! REAL,INTENT(INOUT)         :: Vol(1:dim1,0:N_in,0:N_in,0:N_in) !unchanged volume data, face points will be replaced by vol_in, inner points will be moved by transfinite deformation
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! INTEGER                    :: i,j,k
! REAL                       :: xi0(0:N_in),xiN(0:N_in)
! !REAL                       :: dVol(1:dim1,0:N_in,0:N_in,0:N_in) !Face data, inner points will be overwritten
! !-----------------------------------------------------------------------------------------------------------------------------------
! IF(N_in.EQ.1) THEN
!   Vol=Vol_in
!   RETURN ! no inner nodes
! END IF
! !now N_in>1
! !difference between new (Vol_in) and old (Vold)
! !dVol=Vol_in-Vol !only face data will be used
! Vol=Vol_in

! xiN=0.5*(xi_in+1.)  ![-1,1] ->  0...1
! xi0=1.-xiN          ! 1 ... 0
! !overwrite inner nodes of dVol by transfinite blending of its faces
! DO k=1,N_in-1
!   DO j=1,N_in-1
!     DO i=1,N_in-1
!        Vol(:,i,j,k)=  Vol(:,i,j,0)*xi0(k)+Vol(:,i,j,N_in)*xiN(k)     & !faces
!                       +Vol(:,i,0,k)*xi0(j)+Vol(:,i,N_in,k)*xiN(j)     &
!                       +Vol(:,0,j,k)*xi0(i)+Vol(:,N_in,j,k)*xiN(i)     &
!                       -(  Vol(:,   0,   0,   k)*xi0(i)*xi0(j)          & !-edges
!                         + Vol(:,   0,N_in,   k)*xi0(i)*xiN(j)          &
!                         + Vol(:,   0,   j,   0)*xi0(i)       *xi0(k)   &
!                         + Vol(:,   0,   j,N_in)*xi0(i)       *xiN(k)   &
!                         + Vol(:,N_in,   0,   k)*xiN(i)*xi0(j)          &
!                         + Vol(:,N_in,N_in,   k)*xiN(i)*xiN(j)          &
!                         + Vol(:,N_in,   j,   0)*xiN(i)       *xi0(k)   &
!                         + Vol(:,N_in,   j,N_in)*xiN(i)       *xiN(k)   &
!                         + Vol(:,   i,   0,   0)       *xi0(j)*xi0(k)   &
!                         + Vol(:,   i,   0,N_in)       *xi0(j)*xiN(k)   &
!                         + Vol(:,   i,N_in,   0)       *xiN(j)*xi0(k)   &
!                         + Vol(:,   i,N_in,N_in)       *xiN(j)*xiN(k) ) &
!                       +(  Vol(:,   0,   0,   0)*xi0(i)*xi0(j)*xi0(k)   & ! + corners
!                         + Vol(:,   0,   0,N_in)*xi0(i)*xi0(j)*xiN(k)   &
!                         + Vol(:,   0,N_in,   0)*xi0(i)*xiN(j)*xi0(k)   &
!                         + Vol(:,   0,N_in,N_in)*xi0(i)*xiN(j)*xiN(k)   &
!                         + Vol(:,N_in,   0,   0)*xiN(i)*xi0(j)*xi0(k)   &
!                         + Vol(:,N_in,   0,N_in)*xiN(i)*xi0(j)*xiN(k)   &
!                         + Vol(:,N_in,N_in,   0)*xiN(i)*xiN(j)*xi0(k)   &
!                         + Vol(:,N_in,N_in,N_in)*xiN(i)*xiN(j)*xiN(k) )
!     END DO !i
!   END DO !j
! END DO !k
! !add deformation to volume mapping
! !Vol=Vol+dVol  !vol=vol_in is recovered on the faces, inner nodes will be moved by dVol 
! END SUBROUTINE TransVol


! SUBROUTINE DeformMesh()
! !===================================================================================================================================
! ! Subroutine to read the mesh from a mesh data file
! !===================================================================================================================================
! ! MODULES
! USE MODH_Globals
! USE MODH_Mesh_Vars, ONLY: nTrees,XGeo,Ngeo,Deform
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT VARIABLES
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! OUTPUT VARIABLES
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! INTEGER                        :: i,j,k
! INTEGER                        :: iTree
! REAL                           :: Pi,x(3)
! !-----------------------------------------------------------------------------------------------------------------------------------
! IF(Deform.EQ.0) RETURN
! !deform the mesh
! SELECT CASE(Deform)
! CASE(1) !sinus -1,1 deformation
!   Pi = ACOS(-1.) 
!   DO iTree=1,nTrees
!     DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
!       x(:)=Xgeo(:,i,j,k,iTree)
!       Xgeo(:,i,j,k,iTree) = x+ 0.1*SIN(Pi*x(1))*SIN(Pi*x(2))*SIN(Pi*x(3))
!     END DO; END DO; END DO;
!   END DO
! CASE DEFAULT
!   STOP 'This deform case is not defined'
! END SELECT !Deform
! END SUBROUTINE DeformMesh


SUBROUTINE FinalizeMesh()
!============================================================================================================================
! Deallocate all global interpolation variables.
!============================================================================================================================
! MODULES
! USE MOD_Globals
USE MODH_Mesh_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
!input parameters
!----------------------------------------------------------------------------------------------------------------------------
!output parameters
!----------------------------------------------------------------------------------------------------------------------------
!local variables
INTEGER       :: iTree,iLocSide,iNode
!============================================================================================================================
! Deallocate global variables, needs to go somewhere else later
DO iTree=1,nTrees
  DO iLocSide=1,6
    DEALLOCATE(Trees(iTree)%ep%Side(iLocSide)%sp)
  END DO
  DEALLOCATE(Trees(iTree)%ep)
END DO
DEALLOCATE(Trees)
DO iNode=1,nUniqueNodes
  ADEALLOCATE(UniqueNodes(iNode)%np)
END DO
DEALLOCATE(UniqueNodes)
SDEALLOCATE(XGeo)
SDEALLOCATE(HexMap)
SDEALLOCATE(HexMap_out)
SDEALLOCATE(HexMapInv)
! SDEALLOCATE(Xi_NGeo)
! SDEALLOCATE(XiCL_NGeo_out)
! SDEALLOCATE(wBary_NGeo)
! SDEALLOCATE(Vdm_CL_EQ_out)
SDEALLOCATE(BoundaryName)
SDEALLOCATE(BoundaryType)
END SUBROUTINE FinalizeMesh

END MODULE MODH_Mesh
