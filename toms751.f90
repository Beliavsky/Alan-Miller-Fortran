MODULE TRIPACK
!  ALGORITHM 751, COLLECTED ALGORITHMS FROM ACM.
!  THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!  VOL. 22, NO. 1, March, 1996, P. 1--8.
!  ####### With remark from renka (to appear) 4/dec/1998

! Code converted using TO_F90 by Alan Miller
! Date: 2000-05-25  Time: 18:21:09

IMPLICIT NONE

! Tolerance stored by TRMESH or TRMSHR.
! COMMON/swpcom/swtol
REAL, SAVE, PRIVATE  :: swtol

! REAL :: y
! COMMON/stcom/y
REAL, SAVE, PRIVATE  :: stcom

! Random number seeds
INTEGER, SAVE , PRIVATE :: ix = 1, iy = 2, iz = 3

CONTAINS


SUBROUTINE addcst (ncc, lcc, n, x, y, lwk, list, lptr, lend, ier)

! N.B. Argument IWK has been removed.
 
INTEGER, INTENT(IN)      :: ncc
INTEGER, INTENT(IN)      :: lcc(:)
INTEGER, INTENT(IN)      :: n
REAL, INTENT(IN OUT)     :: x(:)
REAL, INTENT(IN OUT)     :: y(:)
INTEGER, INTENT(IN OUT)  :: lwk
INTEGER, INTENT(IN OUT)  :: list(:)
INTEGER, INTENT(IN OUT)  :: lptr(:)
INTEGER, INTENT(IN OUT)  :: lend(:)
INTEGER, INTENT(OUT)     :: ier

!*******************************************************************

!                                                       From TRIPACK
!                                                    Robert J. Renka
!                                          Dept. of Computer Science
!                                               Univ. of North Texas
!                                                   renka@cs.unt.edu
!                                                           11/12/94

!   This subroutine provides for creation of a constrained Delaunay
! triangulation which, in some sense, covers an arbitrary connected
! region R rather than the convex hull of the nodes.
! This is achieved simply by forcing the presence of certain
! adjacencies (triangulation arcs) corresponding to constraint curves.
! The union of triangles coincides with the convex hull of the nodes,
! but triangles in R can be distinguished from those outside of R.
! The only modification required to generalize the definition of the
! Delaunay triangulation is replacement of property 5 (refer to TRMESH)
! by the following:

!  5')  If a node is contained in the interior of the circum]circle
!       of a triangle, then every interior point of the triangle
!       is separated from the node by a constraint arc.

!   In order to be explicit, we make the following definitions.
! A constraint region is the open interior of a simple closed
! positively oriented polygonal curve defined by an ordered
! sequence of three or more distinct nodes (constraint nodes)
! P(1),P(2),...,P(K), such that P(I) is adjacent to P(I+1) for
! I = 1,...,K with P(K+1) = P(1).
! Thus, the constraint region is on the left (and may have nonfinite
! area) as the sequence of constraint nodes is traversed in the
! specified order.  The constraint regions must not contain nodes
! and must not overlap.  The region R is the convex hull of the nodes
! with constraint regions excluded.

!   Note that the terms boundary node and boundary arc are reserved
! for nodes and arcs on the boundary of the convex hull of the nodes.

!   The algorithm is as follows:  given a triangulation which includes
! one or more sets of constraint nodes, the corresponding adjacencies
! (constraint arcs) are forced to be present (Subroutine EDGE).
! Any additional new arcs required are chosen to be locally optimal
! (satisfy the modified circumcircle property).


! On input:

!       NCC = Number of constraint curves (constraint regions)
!             NCC >= 0.

!       LCC = Array of length NCC (or dummy array of length
!             1 if NCC = 0) containing the index (for X, Y,
!             and LEND) of the first node of constraint I in
!             LCC(I) for I = 1 to NCC.  Thus, constraint I
!             contains K = LCC(I+1) - LCC(I) nodes, K >=
!             3, stored in (X,Y) locations LCC(I), ...,
!             LCC(I+1)-1, where LCC(NCC+1) = N+1.

!       N = Number of nodes in the triangulation, including
!           constraint nodes.  N >= 3.

!       X,Y = Arrays of length N containing the coordinates of the nodes with
!             non-constraint nodes in the first LCC(1)-1 locations, followed
!             by NCC sequences of constraint nodes.  Only one of these
!             sequences may be specified in clockwise order to represent an
!             exterior constraint curve (a constraint region with nonfinite
!             area).

! The above parameters are not altered by this routine.

!       LWK = Length of IWK.  This must be at least 2*NI where NI is the
!             maximum number of arcs which intersect a constraint arc to be
!             added.  NI is bounded by N-3.

!       IWK = Integer work array of length LWK (used by
!             Subroutine EDGE to add constraint arcs).

!       LIST,LPTR,LEND = Data structure defining the triangulation.
!                        Refer to Subroutine TRMESH.

! On output:

!       LWK = Required length of IWK unless IER = 1 or IER = 3.
!             In the case of IER = 1, LWK is not altered
!             from its input value.

!       IWK = Array containing the endpoint indexes of the new arcs
!             which were swapped in by the last call to Subroutine EDGE.

!       LIST,LPTR,LEND = Triangulation data structure with all constraint
!                        arcs present unless  IER .NE. 0.
!                        These arrays are not altered if IER = 1.

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if NCC, N, or an LCC entry is outside
!                     its valid range, or LWK < 0 on input.
!             IER = 2 if more space is required in IWK.
!             IER = 3 if the triangulation data structure is invalid,
!                     or failure (in EDGE or OPTIM) was caused by
!                     collinear nodes on the convex hull boundary.
!                     An error message is written to logical unit 6
!                     in this case.
!             IER = 4 if intersecting constraint arcs were encountered.
!             IER = 5 if a constraint region contains a node.

! Modules required by ADDCST:  EDGE, LEFT, LSTPTR, OPTIM, SWAP, SWPTST

! Intrinsic functions called by ADDCST:  ABS, MAX

!***********************************************************

INTEGER :: i, ifrst, ilast, k, kbak, kfor, kn, lccip1,  &
           lp, lpb, lpf, lpl, lw, lwd2, n1, n2

lwd2 = lwk/2

! Test for errors in input parameters.

ier = 1
IF (ncc < 0 .OR. lwk < 0) RETURN
IF (ncc == 0) THEN
  IF (n < 3) RETURN
  lwk = 0
  GO TO 9
ELSE
  lccip1 = n+1
  DO  i = ncc,1,-1
    IF (lccip1 - lcc(i) < 3) RETURN
    lccip1 = lcc(i)
  END DO
  IF (lccip1 < 1) RETURN
END IF

! Force the presence of constraint arcs.  The outer loop is
!   on constraints in reverse order.  IFRST and ILAST are
!   the first and last nodes of constraint I.

lwk = 0
ifrst = n+1
DO  i = ncc,1,-1
  ilast = ifrst - 1
  ifrst = lcc(i)
  
!   Inner loop on constraint arcs N1-N2 in constraint I.
  
  n1 = ilast
  DO  n2 = ifrst,ilast
    lw = lwd2
    CALL edge (n1,n2,x,y, lw,list,lptr,lend, ier)
    lwk = MAX(lwk,2*lw)
    IF (ier == 4) ier = 3
    IF (ier /= 0) RETURN
    n1 = n2
  END DO
END DO

! Test for errors.  The outer loop is on constraint I with first and
!   last nodes IFRST and ILAST, and the inner loop is on constraint nodes K
!   with (KBAK,K,KFOR) a subsequence of constraint I.

ier = 4
ifrst = n+1
DO  i = ncc,1,-1
  ilast = ifrst - 1
  ifrst = lcc(i)
  kbak = ilast
  DO  k = ifrst,ilast
    kfor = k + 1
    IF (k == ilast) kfor = ifrst
    
!   Find the LIST pointers LPF and LPB of KFOR and KBAK as neighbors of K.
    
    lpf = 0
    lpb = 0
    lpl = lend(k)
    lp = lpl

    DO
      lp = lptr(lp)
      kn = ABS(list(lp))
      IF (kn == kfor) lpf = lp
      IF (kn == kbak) lpb = lp
      IF (lp == lpl) EXIT
    END DO
    
!   A pair of intersecting constraint arcs was encountered
!     if and only if a constraint arc is missing (introduction
!     of the second caused the first to be swapped out).
    
    IF (lpf == 0 .OR. lpb == 0) RETURN
    
!   Loop on neighbors KN of node K which follow KFOR and precede KBAK.
!     The constraint region contains no nodes
!     if and only if all such nodes KN are in constraint I.
    
    lp = lpf
    DO
      lp = lptr(lp)
      IF (lp == lpb) EXIT
      kn = ABS(list(lp))
      IF (kn < ifrst .OR. kn > ilast) GO TO 10
    END DO
    
!   Bottom of loop.
    
    kbak = k
  END DO
END DO

! No errors encountered.

9 ier = 0
RETURN

! A constraint region contains a node.

10 ier = 5
RETURN
END SUBROUTINE addcst



SUBROUTINE addnod (k, xk, yk, ist, ncc, lcc, n, x, y, list, lptr, lend,  &
                   lnew, ier)

INTEGER, INTENT(IN)      :: k
REAL, INTENT(IN)         :: xk
REAL, INTENT(IN)         :: yk
INTEGER, INTENT(IN)      :: ist
INTEGER, INTENT(IN)      :: ncc
INTEGER, INTENT(IN OUT)  :: lcc(:)
INTEGER, INTENT(IN OUT)  :: n
REAL, INTENT(IN OUT)     :: x(:)
REAL, INTENT(IN OUT)     :: y(:)
INTEGER, INTENT(IN OUT)  :: list(:)
INTEGER, INTENT(IN OUT)  :: lptr(:)
INTEGER, INTENT(IN OUT)  :: lend(:)
INTEGER, INTENT(IN OUT)  :: lnew
INTEGER, INTENT(OUT)     :: ier


!*******************************************************************

!                                                       From TRIPACK
!                                                    Robert J. Renka
!                                          Dept. of Computer Science
!                                               Univ. of North Texas
!                                                   renka@cs.unt.edu
!                                                           06/27/98

!   Given a triangulation of N nodes in the plane created by
! Subroutine TRMESH or TRMSHR, this subroutine updates the data
! structure with the addition of a new node in position K.
! If node K is inserted into X and Y (K <= N) rather than appended
! (K = N+1), then a corresponding insertion must be performed in
! any additional arrays associated with the nodes.
! For example, an array of data values Z must be shifted down
! to open up position K for the new value:
! set Z(I+1) to Z(I) for I = N,N-1,...,K.  For optimal efficiency,
! new nodes should be appended whenever possible.
! Insertion is necessary, however, to add a non-constraint node
! when constraints are present (refer to Subroutine ADDCST).

!   Note that a constraint node cannot be added by this routine.
! In order to insert a constraint node, it is necessary to add the node with
! no constraints present (call this routine with NCC = 0), update LCC by
! incrementing the appropriate entries, and then create (or restore) the
! constraints by a call to ADDCST.

!   The algorithm consists of the following steps:  node K is located
! relative to the triangulation (TRFIND), its index is added to the data
! structure (INTADD or BDYADD), and a sequence of swaps (SWPTST and SWAP) are
! applied to the arcs opposite K so that all arcs incident on node K and
! opposite node K (excluding constraint arcs) are locally optimal (satisfy
! the circumcircle test).  Thus, if a (constrained) Delaunay triangulation is
! input, a (constrained) Delaunay triangulation will result.
! All indexes are incremented as necessary for an insertion.

! On input:

!       K = Nodal index (index for X, Y, and LEND) of the
!           new node to be added.  1 <= K <= LCC(1).
!           (K <= N+1 if NCC=0).

!       XK,YK = Cartesian coordinates of the new node (to be stored in X(K)
!               and Y(K)).  The node must not lie in a constraint region.

!       IST = Index of a node at which TRFIND begins the search.
!             Search time depends on the proximity of this node to node K.
!             1 <= IST <= N.

!       NCC = Number of constraint curves.  NCC >= 0.

! The above parameters are not altered by this routine.

!       LCC = List of constraint curve starting indexes (or dummy array of
!             length 1 if NCC = 0).  Refer to Subroutine ADDCST.

!       N = Number of nodes in the triangulation before K is added.  N >= 3.
!           Note that N will be incremented following the addition of node K.

!       X,Y = Arrays of length at least N+1 containing the Cartesian
!             coordinates of the nodes in the first N positions with non-
!             constraint nodes in the first LCC(1)-1 locations if NCC > 0.

!       LIST,LPTR,LEND,LNEW = Data structure associated with the triangulation
!                             of nodes 1 to N.  The arrays must have
!                             sufficient length for N+1 nodes.
!                             Refer to TRMESH.

! On output:

!       LCC = List of constraint curve starting indexes incremented by 1
!             to reflect the insertion of K unless NCC = 0 or
!             (IER .NE. 0 and IER .NE. -4).

!       N = Number of nodes in the triangulation including K
!           unless IER .NE. 0 and IER .NE. -4.  Note that
!           all comments refer to the input value of N.

!       X,Y = Arrays updated with the insertion of XK and YK in the K-th
!             positions (node I+1 was node I be fore the insertion for
!             I = K to N if K <= N) unless IER .NE. 0 and IER .NE. -4.

!       LIST,LPTR,LEND,LNEW = Data structure updated with the addition of
!                             node K unless IER .NE. 0 and IER .NE. -4.

!       IER = Error indicator:
!             IER =  0 if no errors were encountered.
!             IER = -1 if K, IST, NCC, N, or an LCC entry is
!                      outside its valid range on input.
!             IER = -2 if all nodes (including K) are collinear.
!             IER =  L if nodes L and K coincide for some L.
!             IER = -3 if K lies in a constraint region.
!             IER = -4 if an error flag is returned by SWAP implying that the
!                      triangulation (geometry) was bad on input.

!             The errors conditions are tested in the order specified.

! Modules required by ADDNOD:  BDYADD, CRTRI, INDXCC, INSERT, INTADD, JRAND,
!                                LEFT, LSTPTR, SWAP, SWPTST, TRFIND

! Intrinsic function called by ADDNOD:  ABS

!***********************************************************

INTEGER  :: i, i1, i2, i3, ibk, io1, io2, in1, kk, l, lccip1, lp, lpf,  &
            lpo1, nm1

kk = k

! Test for an invalid input parameter.

IF (kk < 1 .OR. ist < 1 .OR. ist > n .OR.  ncc < 0 .OR. n < 3) GO TO 7
lccip1 = n+1
DO  i = ncc,1,-1
  IF (lccip1-lcc(i) < 3) GO TO 7
  lccip1 = lcc(i)
END DO
IF (kk > lccip1) GO TO 7

! Find a triangle (I1,I2,I3) containing K or the rightmost (I1)
!   and leftmost (I2) visible boundary nodes as viewed from node K.

CALL trfind (ist, xk, yk, n, x, y, list, lptr, lend, i1, i2, i3)

! Test for collinear nodes, duplicate nodes, and K lying in
!   a constraint region.

IF (i1 == 0) GO TO 8
IF (i3 /= 0) THEN
  l = i1
  IF (xk == x(l) .AND. yk == y(l)) GO TO 9
  l = i2
  IF (xk == x(l) .AND. yk == y(l)) GO TO 9
  l = i3
  IF (xk == x(l) .AND. yk == y(l)) GO TO 9
  IF (ncc > 0 .AND. crtri(ncc, lcc, i1, i2, i3) ) GO TO 10
ELSE
  
!   K is outside the convex hull of the nodes and lies in a
!     constraint region iff an exterior constraint curve is present.
  
  IF (ncc > 0 .AND. indxcc(ncc,lcc,n,list,lend)  /= 0) GO TO 10
END IF

! No errors encountered.

ier = 0
nm1 = n
n = n + 1
IF (kk < n) THEN
  
! Open a slot for K in X, Y, and LEND, and increment all
!   nodal indexes which are greater than or equal to K.
!   Note that LIST, LPTR, and LNEW are not yet updated with
!   either the neighbors of K or the edges terminating on K.
  
  DO  ibk = nm1,kk,-1
    x(ibk+1) = x(ibk)
    y(ibk+1) = y(ibk)
    lend(ibk+1) = lend(ibk)
  END DO
  DO  i = 1,ncc
    lcc(i) = lcc(i) + 1
  END DO
  l = lnew - 1
  DO  i = 1,l
    IF (list(i) >= kk) list(i) = list(i) + 1
    IF (list(i) <= -kk) list(i) = list(i) - 1
  END DO
  IF (i1 >= kk) i1 = i1 + 1
  IF (i2 >= kk) i2 = i2 + 1
  IF (i3 >= kk) i3 = i3 + 1
END IF

! Insert K into X and Y, and update LIST, LPTR, LEND, and
!   LNEW with the arcs containing node K.

x(kk) = xk
y(kk) = yk
IF (i3 == 0) THEN
  CALL bdyadd (kk, i1, i2, list, lptr, lend, lnew )
ELSE
  CALL intadd (kk, i1, i2, i3, list, lptr, lend, lnew )
END IF

! Initialize variables for optimization of the triangulation.

lp = lend(kk)
lpf = lptr(lp)
io2 = list(lpf)
lpo1 = lptr(lpf)
io1 = ABS(list(lpo1))

! Begin loop:  find the node opposite K.

5 lp = lstptr(lend(io1),io2,list,lptr)
IF (list(lp) < 0) GO TO 6
lp = lptr(lp)
in1 = ABS(list(lp))
IF ( crtri(ncc,lcc,io1,io2,in1) ) GO TO 6

! Swap test: if a swap occurs, two new arcs are opposite K and must be tested.

IF ( .NOT. swptst(in1,kk,io1,io2,x,y) ) GO TO 6
CALL swap (in1,kk,io1,io2, list,lptr,lend, lpo1)
IF (lpo1 == 0) GO TO 11
io1 = in1
GO TO 5

! No swap occurred.  Test for termination and reset IO2 and IO1.

6 IF (lpo1 == lpf .OR. list(lpo1) < 0) RETURN
io2 = io1
lpo1 = lptr(lpo1)
io1 = ABS(list(lpo1))
GO TO 5

! A parameter is outside its valid range on input.

7 ier = -1
RETURN

! All nodes are collinear.

8 ier = -2
RETURN

! Nodes L and K coincide.

9 ier = l
RETURN

! Node K lies in a constraint region.

10 ier = -3
RETURN

! Zero pointer returned by SWAP.

11 ier = -4
RETURN
END SUBROUTINE addnod



FUNCTION areap (x, y, nb, nodes) RESULT(fn_val)

REAL, INTENT(IN)     :: x(:)
REAL, INTENT(IN)     :: y(:)
INTEGER, INTENT(IN)  :: nb
INTEGER, INTENT(IN)  :: nodes(:)
REAL                 :: fn_val

!*******************************************************************

!                                                       From TRIPACK
!                                                    Robert J. Renka
!                                          Dept. of Computer Science
!                                               Univ. of North Texas
!                                                   renka@cs.unt.edu
!                                                           09/21/90

!   Given a sequence of NB points in the plane, this function computes
! the signed area bounded by the closed polygonal curve which passes
! through the points in the specified order.
! Each simple closed curve is positively oriented (bounds positive area)
! if and only if the points are specified in counterclockwise order.
! The last point of the curve is taken to be the first point specified,
! and this point should therefore not be specified twice.

!   The area of a triangulation may be computed by calling
! AREAP with values of NB and NODES determined by Subroutine BNODES.


! On input:

!       X,Y = Arrays of length N containing the Cartesian coordinates of a
!             set of points in the plane for some N >= NB.

!       NB = Length of NODES.

!       NODES = Array of length NB containing the ordered sequence of nodal
!               indexes (in the range 1 to N) which define the polygonal curve.

! Input parameters are not altered by this function.

! On output:

!       AREAP = Signed area bounded by the polygonal curve, or zero if NB < 3.

! Modules required by AREAP:  None

!***********************************************************

INTEGER :: i, nd1, nd2, nnb
REAL    :: a

! Local parameters:

! A =       Partial sum of signed (and doubled) trapezoid areas
! I =       DO-loop and NODES index
! ND1,ND2 = Elements of NODES
! NNB =     Local copy of NB

nnb = nb
a = 0.
IF (nnb < 3) GO TO 2
nd2 = nodes(nnb)

! Loop on line segments NODES(I-1) -> NODES(I), where
!   NODES(0) = NODES(NB), adding twice the signed trapezoid
!   areas (integrals of the linear interpolants) to A.

DO  i = 1,nnb
  nd1 = nd2
  nd2 = nodes(i)
  a = a + (x(nd2)-x(nd1))*(y(nd1)+y(nd2))
END DO

! A contains twice the negative signed area of the region.

2 fn_val = -a/2.
RETURN
END FUNCTION areap




SUBROUTINE bdyadd (kk, i1, i2, list, lptr, lend, lnew )

INTEGER, INTENT(IN)      :: kk
INTEGER, INTENT(IN)      :: i1
INTEGER, INTENT(IN)      :: i2
INTEGER, INTENT(OUT)     :: list(:)
INTEGER, INTENT(IN OUT)  :: lptr(:)
INTEGER, INTENT(IN OUT)  :: lend(:)
INTEGER, INTENT(IN OUT)  :: lnew

!*******************************************************************

!                                                       From TRIPACK
!                                                    Robert J. Renka
!                                          Dept. of Computer Science
!                                               Univ. of North Texas
!                                                   renka@cs.unt.edu
!                                                           02/22/91

!   This subroutine adds a boundary node to a triangulation of a set
! of points in the plane.  The data structure is updated with the
! insertion of node KK, but no optimization is performed.


! On input:

!       KK = Index of a node to be connected to the sequence
!            of all visible boundary nodes.  KK >= 1 and
!            KK must not be equal to I1 or I2.

!       I1 = First (rightmost as viewed from KK) boundary
!            node in the triangulation which is visible from
!            node KK (the line segment KK-I1 intersects no arcs.

!       I2 = Last (leftmost) boundary node which is visible from node KK.
!            I1 and I2 may be determined by Subroutine TRFIND.

! The above parameters are not altered by this routine.

!       LIST,LPTR,LEND,LNEW = Triangulation data structure created by TRMESH
!                             or TRMSHR.  Nodes I1 and I2 must be included
!                             in the triangulation.

! On output:

!       LIST,LPTR,LEND,LNEW = Data structure updated with the addition of
!                             node KK.  Node KK is connected to I1, I2, and
!                             all boundary nodes in between.

! Module required by BDYADD:  INSERT

!***********************************************************

INTEGER :: k, lp, lsav, n1, n2, next, nsav

k = kk
n1 = i1
n2 = i2

! Add K as the last neighbor of N1.

lp = lend(n1)
lsav = lptr(lp)
lptr(lp) = lnew
list(lnew) = -k
lptr(lnew) = lsav
lend(n1) = lnew
lnew = lnew + 1
next = -list(lp)
list(lp) = next
nsav = next

! Loop on the remaining boundary nodes between N1 and N2,
!   adding K as the first neighbor.

1 lp = lend(next)
CALL insert (k,lp,list,lptr,lnew)
IF (next == n2) GO TO 2
next = -list(lp)
list(lp) = next
GO TO 1

! Add the boundary nodes between N1 and N2 as neighbors of node K.

2 lsav = lnew
list(lnew) = n1
lptr(lnew) = lnew + 1
lnew = lnew + 1
next = nsav

3 IF (next == n2) GO TO 4
list(lnew) = next
lptr(lnew) = lnew + 1
lnew = lnew + 1
lp = lend(next)
next = list(lp)
GO TO 3

4 list(lnew) = -n2
lptr(lnew) = lsav
lend(k) = lnew
lnew = lnew + 1
RETURN
END SUBROUTINE bdyadd



SUBROUTINE bnodes (n, list, lptr, lend, nodes, nb, na, nt)

INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: list(:)
INTEGER, INTENT(IN)   :: lptr(:)
INTEGER, INTENT(IN)   :: lend(:)
INTEGER, INTENT(OUT)  :: nodes(:)
INTEGER, INTENT(OUT)  :: nb
INTEGER, INTENT(OUT)  :: na
INTEGER, INTENT(OUT)  :: nt

!*******************************************************************

!                                                       From TRIPACK
!                                                    Robert J. Renka
!                                          Dept. of Computer Science
!                                               Univ. of North Texas
!                                                   renka@cs.unt.edu
!                                                           09/01/88

!   Given a triangulation of N points in the plane, this subroutine
! returns an array containing the indexes, in counterclockwise order,
! of the nodes on the boundary of the convex hull of the set of points.


! On input:

!       N = Number of nodes in the triangulation.  N >= 3.

!       LIST,LPTR,LEND = Data structure defining the triangulation.
!                        Refer to Subroutine TRMESH.

! The above parameters are not altered by this routine.

!       NODES = Integer array of length at least NB
!               (NB <= N).

! On output:

!       NODES = Ordered sequence of boundary node indexes in the range 1 to N.

!       NB = Number of boundary nodes.

!       NA,NT = Number of arcs and triangles, respectively,
!               in the triangulation.

! Modules required by BNODES:  None

!***********************************************************

INTEGER :: k, lp, n0, nst

! Set NST to the first boundary node encountered.

nst = 1
1 lp = lend(nst)
IF (list(lp) < 0) GO TO 2
nst = nst + 1
GO TO 1

! Initialization.

2 nodes(1) = nst
k = 1
n0 = nst

! Traverse the boundary in counterclockwise order.

3 lp = lend(n0)
lp = lptr(lp)
n0 = list(lp)
IF (n0 == nst) GO TO 4
k = k + 1
nodes(k) = n0
GO TO 3

! Termination.

4 nb = k
nt = 2*n - nb - 2
na = nt + n - 1
RETURN
END SUBROUTINE bnodes




SUBROUTINE circum (x1, y1, x2, y2, x3, y3, ratio, xc, yc, cr, sa, ar)

REAL, INTENT(IN)     :: x1
REAL, INTENT(IN)     :: y1
REAL, INTENT(IN)     :: x2
REAL, INTENT(IN)     :: y2
REAL, INTENT(IN)     :: x3
REAL, INTENT(IN)     :: y3
LOGICAL, INTENT(IN)  :: ratio
REAL, INTENT(OUT)    :: xc
REAL, INTENT(OUT)    :: yc
REAL, INTENT(OUT)    :: cr
REAL, INTENT(OUT)    :: sa
REAL, INTENT(OUT)    :: ar

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   12/10/96

!   Given three vertices defining a triangle, this subroutine returns
! the circumcenter, circumradius, signed triangle area, and,
! optionally, the aspect ratio of the triangle.


! On input:

!       X1,...,Y3 = Cartesian coordinates of the vertices.

!       RATIO = Logical variable with value TRUE if and only
!               if the aspect ratio is to be computed.

! Input parameters are not altered by this routine.

! On output:

!       XC,YC = Cartesian coordinates of the circumcenter
!               (center of the circle defined by the three points)
!               unless SA = 0, in which XC and YC are not altered.

!       CR = Circumradius (radius of the circle defined by
!            the three points) unless SA = 0 (infinite
!            radius), in which case CR is not altered.

!       SA = Signed triangle area with positive value if and only if
!            the vertices are specified in counterclockwise order:
!            (X3,Y3) is strictly to the left of the directed line
!            from (X1,Y1) toward (X2,Y2).

!       AR = Aspect ratio r/CR, where r is the radius of the inscribed
!            circle, unless RATIO = FALSE, in which case AR is not
!            altered.  AR is in the range 0 to .5, with value 0 iff
!            SA = 0 and value .5 iff the vertices define an equilateral
!            triangle.

! Modules required by CIRCUM:  None

! Intrinsic functions called by CIRCUM:  ABS, SQRT

!***********************************************************

INTEGER :: i
REAL    :: ds(3), fx, fy, u(3), v(3)

! Set U(K) and V(K) to the x and y components, respectively,
!   of the directed edge opposite vertex K.

u(1) = x3 - x2
u(2) = x1 - x3
u(3) = x2 - x1
v(1) = y3 - y2
v(2) = y1 - y3
v(3) = y2 - y1

! Set SA to the signed triangle area.

sa = (u(1)*v(2) - u(2)*v(1))/2.
IF (sa == 0.) THEN
  IF (ratio) ar = 0.
  RETURN
END IF

! Set DS(K) to the squared distance from the origin to vertex K.

ds(1) = x1*x1 + y1*y1
ds(2) = x2*x2 + y2*y2
ds(3) = x3*x3 + y3*y3

! Compute factors of XC and YC.

fx = 0.
fy = 0.
DO  i = 1,3
  fx = fx - ds(i)*v(i)
  fy = fy + ds(i)*u(i)
END DO
xc = fx/(4.*sa)
yc = fy/(4.*sa)
cr = SQRT( (xc-x1)**2 + (yc-y1)**2 )
IF (.NOT. ratio) RETURN

! Compute the squared edge lengths and aspect ratio.

DO  i = 1,3
  ds(i) = u(i)*u(i) + v(i)*v(i)
END DO
ar = 2.*ABS(sa)/ ( (SQRT(ds(1)) + SQRT(ds(2)) + SQRT(ds(3)))*cr )
RETURN
END SUBROUTINE circum



FUNCTION crtri (ncc, lcc, i1, i2, i3) RESULT(fn_val)

INTEGER, INTENT(IN)  :: ncc
INTEGER, INTENT(IN)  :: lcc(:)
INTEGER, INTENT(IN)  :: i1
INTEGER, INTENT(IN)  :: i2
INTEGER, INTENT(IN)  :: i3
LOGICAL              :: fn_val

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   08/14/91

!   This function returns TRUE if and only if triangle (I1,
! I2,I3) lies in a constraint region.


! On input:

!       NCC,LCC = Constraint data structure.  Refer to Subroutine ADDCST.

!       I1,I2,I3 = Nodal indexes of the counterclockwise-
!                  ordered vertices of a triangle.

! Input parameters are altered by this function.

!       CRTRI = TRUE iff (I1,I2,I3) is a constraint region triangle.

! Note that input parameters are not tested for validity.

! Modules required by CRTRI:  None

! Intrinsic functions called by CRTRI:  MAX, MIN

!***********************************************************

INTEGER :: i, imax, imin

imax = MAX(i1, i2, i3)

!   Find the index I of the constraint containing IMAX.

i = ncc + 1
1 i = i - 1
IF (i <= 0) GO TO 2
IF (imax < lcc(i)) GO TO 1
imin = MIN(i1, i2, i3)

! P lies in a constraint region iff I1, I2, and I3 are nodes of the same
!   constraint (IMIN >= LCC(I)), and (IMIN,IMAX) is (I1,I3), (I2,I1),
!   or (I3,I2).

fn_val = imin >= lcc(i) .AND. ((imin == i1 .AND. imax == i3) .OR. &
         (imin == i2 .AND. imax == i1) .OR. (imin == i3 .AND. imax == i2))
RETURN

! NCC <= 0 or all vertices are non-constraint nodes.

2 fn_val = .false.
RETURN
END FUNCTION crtri



SUBROUTINE delarc (n, io1, io2, list, lptr, lend, lnew, ier)

INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN)      :: io1
INTEGER, INTENT(IN)      :: io2
INTEGER, INTENT(IN OUT)  :: list(:)
INTEGER, INTENT(IN OUT)  :: lptr(:)
INTEGER, INTENT(IN OUT)  :: lend(:)
INTEGER, INTENT(IN OUT)  :: lnew
INTEGER, INTENT(OUT)     :: ier

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   11/12/94

!   This subroutine deletes a boundary arc from a triangulation.
! It may be used to remove a null triangle from the convex hull boundary.
! Note, however, that if the union of triangles is rendered nonconvex,
! Subroutines DELNOD, EDGE, and TRFIND may fail.
! Thus, Subroutines ADDCST, ADDNOD, DELNOD, EDGE, and NEARND should not be
! called following an arc deletion.

! On input:

!       N = Number of nodes in the triangulation.  N >= 4.

!       IO1,IO2 = Indexes (in the range 1 to N) of a pair of adjacent
!                 boundary nodes defining the arc to be removed.

! The above parameters are not altered by this routine.

!       LIST,LPTR,LEND,LNEW = Triangulation data structure created by
!                             TRMESH or TRMSHR.

! On output:

!       LIST,LPTR,LEND,LNEW = Data structure updated with the removal of
!                             arc IO1-IO2 unless IER > 0.

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, IO1, or IO2 is outside its valid range,
!                     or IO1 = IO2.
!             IER = 2 if IO1-IO2 is not a boundary arc.
!             IER = 3 if the node opposite IO1-IO2 is already a boundary node,
!                     and thus IO1 or IO2 has only two neighbors or a
!                     deletion would result in two triangulations sharing
!                     a single node.
!             IER = 4 if one of the nodes is a neighbor of the other,
!                     but not vice versa, implying an invalid triangulation
!                     data structure.

! Modules required by DELARC:  DELNB, LSTPTR

! Intrinsic function called by DELARC:  ABS

!***********************************************************

INTEGER :: lp, lph, lpl, n1, n2, n3

n1 = io1
n2 = io2

! Test for errors, and set N1->N2 to the directed boundary edge
!   associated with IO1-IO2:  (N1,N2,N3) is a triangle for some N3.

IF (n < 4 .OR. n1 < 1 .OR. n1 > n .OR. n2 < 1 .OR. n2 > n .OR. n1 == n2) THEN
  ier = 1
  RETURN
END IF

lpl = lend(n2)
IF (-list(lpl) /= n1) THEN
  n1 = n2
  n2 = io1
  lpl = lend(n2)
  IF (-list(lpl) /= n1) THEN
    ier = 2
    RETURN
  END IF
END IF

! Set N3 to the node opposite N1->N2 (the second neighbor
!   of N1), and test for error 3 (N3 already a boundary node).

lpl = lend(n1)
lp = lptr(lpl)
lp = lptr(lp)
n3 = ABS(list(lp))
lpl = lend(n3)
IF (list(lpl) <= 0) THEN
  ier = 3
  RETURN
END IF

! Delete N2 as a neighbor of N1, making N3 the first neighbor, and test for
!   error 4 (N2 not a neighbor of N1).
!   Note that previously computed pointers may no longer be valid following
!   the call to DELNB.

CALL delnb (n1,n2,n, list,lptr,lend,lnew, lph)
IF (lph < 0) THEN
  ier = 4
  RETURN
END IF

! Delete N1 as a neighbor of N2, making N3 the new last neighbor.

CALL delnb (n2,n1,n, list,lptr,lend,lnew, lph)

! Make N3 a boundary node with first neighbor N2 and last neighbor N1.

lp = lstptr(lend(n3), n1, list, lptr)
lend(n3) = lp
list(lp) = -n1

! No errors encountered.

ier = 0
RETURN
END SUBROUTINE delarc



SUBROUTINE delnb (n0, nb, n, list, lptr, lend, lnew, lph)

INTEGER, INTENT(IN)      :: n0
INTEGER, INTENT(IN)      :: nb
INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN OUT)  :: list(:)
INTEGER, INTENT(IN OUT)  :: lptr(:)
INTEGER, INTENT(IN OUT)  :: lend(:)
INTEGER, INTENT(IN OUT)  :: lnew
INTEGER, INTENT(OUT)     :: lph

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/30/98

!   This subroutine deletes a neighbor NB from the adjacency list of node N0
! (but N0 is not deleted from the adjacency list of NB) and, if NB is a
! boundary node, makes N0 a boundary node.  For pointer (LIST index) LPH to NB
! as a neighbor of N0, the empty LIST,LPTR location LPH is filled in with the
! values at LNEW-1, pointer LNEW-1 (in LPTR and possibly in LEND) is changed
! to LPH, and LNEW is decremented.  This requires a search of LEND and LPTR
! entailing an expected operation count of O(N).

! On input:

!       N0,NB = Indexes, in the range 1 to N, of a pair of
!               nodes such that NB is a neighbor of N0.
!               (N0 need not be a neighbor of NB.)

!       N = Number of nodes in the triangulation.  N >= 3.

! The above parameters are not altered by this routine.

!       LIST,LPTR,LEND,LNEW = Data structure defining the triangulation.

! On output:

!       LIST,LPTR,LEND,LNEW = Data structure updated with the removal of
!                             NB from the adjacency list of N0 unless LPH < 0.

!       LPH = List pointer to the hole (NB as a neighbor of
!             N0) filled in by the values at LNEW-1 or error indicator:
!             LPH > 0 if no errors were encountered.
!             LPH = -1 if N0, NB, or N is outside its valid range.
!             LPH = -2 if NB is not a neighbor of N0.

! Modules required by DELNB:  None

! Intrinsic function called by DELNB:  ABS

!***********************************************************

INTEGER :: i, lnw, lp, lpb, lpl, lpp, nn

! Local parameters:

! I =   DO-loop index
! LNW = LNEW-1 (output value of LNEW)
! LP =  LIST pointer of the last neighbor of NB
! LPB = Pointer to NB as a neighbor of N0
! LPL = Pointer to the last neighbor of N0
! LPP = Pointer to the neighbor of N0 that precedes NB
! NN =  Local copy of N

nn = n

! Test for error 1.

IF (n0 < 1 .OR. n0 > nn .OR. nb < 1  .OR. nb > nn .OR. nn < 3) THEN
  lph = -1
  RETURN
END IF

!   Find pointers to neighbors of N0:

!     LPL points to the last neighbor,
!     LPP points to the neighbor NP preceding NB, and
!     LPB points to NB.

lpl = lend(n0)
lpp = lpl
lpb = lptr(lpp)
1 IF (list(lpb) == nb) GO TO 2
lpp = lpb
lpb = lptr(lpp)
IF (lpb /= lpl) GO TO 1

!   Test for error 2 (NB not found).

IF (ABS(list(lpb)) /= nb) THEN
  lph = -2
  RETURN
END IF

!   NB is the last neighbor of N0.  Make NP the new last neighbor and,
!     if NB is a boundary node, then make N0 a boundary node.

lend(n0) = lpp
lp = lend(nb)
IF (list(lp) < 0) list(lpp) = -list(lpp)
GO TO 3

!   NB is not the last neighbor of N0.  If NB is a boundary node and
!     N0 is not, then make N0 a boundary node with last neighbor NP.

2 lp = lend(nb)
IF (list(lp) < 0 .AND. list(lpl) > 0) THEN
  lend(n0) = lpp
  list(lpp) = -list(lpp)
END IF

!   Update LPTR so that the neighbor following NB now follows NP,
!     and fill in the hole at location LPB.

3 lptr(lpp) = lptr(lpb)
lnw = lnew-1
list(lpb) = list(lnw)
lptr(lpb) = lptr(lnw)
DO  i = nn,1,-1
  IF (lend(i) == lnw) THEN
    lend(i) = lpb
    GO TO 5
  END IF
END DO

5 DO  i = 1,lnw-1
  IF (lptr(i) == lnw) THEN
    lptr(i) = lpb
  END IF
END DO

! No errors encountered.

lnew = lnw
lph = lpb
RETURN
END SUBROUTINE delnb



SUBROUTINE delnod (k, ncc, lcc, n, x, y, list, lptr, lend, lnew, lwk, iwk, ier)

INTEGER, INTENT(IN)      :: k
INTEGER, INTENT(IN)      :: ncc
INTEGER, INTENT(IN OUT)  :: lcc(:)
INTEGER, INTENT(IN OUT)  :: n
REAL, INTENT(IN OUT)     :: x(:)
REAL, INTENT(IN OUT)     :: y(:)
INTEGER, INTENT(IN OUT)  :: list(:)
INTEGER, INTENT(IN OUT)  :: lptr(:)
INTEGER, INTENT(IN OUT)  :: lend(:)
INTEGER, INTENT(IN OUT)  :: lnew
INTEGER, INTENT(IN OUT)  :: lwk
INTEGER, INTENT(OUT)     :: iwk(:,:)
INTEGER, INTENT(OUT)     :: ier

!*******************************************************************

!                                                       From TRIPACK
!                                                    Robert J. Renka
!                                          Dept. of Computer Science
!                                               Univ. of North Texas
!                                                   renka@cs.unt.edu
!                                                           06/28/98

!   This subroutine deletes node K (along with all arcs incident on
! node K) from a triangulation of N nodes in the plane, and inserts arcs
! as necessary to produce a triangulation of the remaining N-1 nodes.
! If a Delaunay triangulation is input, a Delaunay triangulation will result,
! and thus, DELNOD reverses the effect of a call to Subroutine ADDNOD.

!   Note that a constraint node cannot be deleted by this routine.
! In order to delete a constraint node, it is necessary to call this
! routine with NCC = 0, decrement the appropriate LCC entries (LCC(I)
! such that LCC(I) > K), and then create (or restore) the constraints
! by a call to Subroutine ADDCST.


! On input:

!       K = Index (for X and Y) of the node to be deleted.
!           1 <= K < LCC(1).  (K <= N if NCC=0).

!       NCC = Number of constraint curves.  NCC >= 0.

! The above parameters are not altered by this routine.

!       LCC = List of constraint curve starting indexes (or dummy array of
!             length 1 if NCC = 0).  Refer to Subroutine ADDCST.

!       N = Number of nodes in the triangulation on input.
!           N >= 4.  Note that N will be decremented following the deletion.

!       X,Y = Arrays of length N containing the coordinates of the nodes with
!             non-constraint nodes in the first LCC(1)-1 locations if NCC > 0.

!       LIST,LPTR,LEND,LNEW = Data structure defining the triangulation.
!                             Refer to Subroutine TRMESH.

!       LWK = Number of columns reserved for IWK.  LWK must be at least NNB-3,
!             where NNB is the number of neighbors of node K, including an
!             extra pseudo-node if K is a boundary node.

!       IWK = Integer work array dimensioned 2 by LWK (or array of
!             length >= 2*LWK).

! On output:

!       LCC = List of constraint curve starting indexes decremented by 1 to
!             reflect the deletion of K unless NCC = 0 or 1 <= IER <= 4.

!       N = New number of nodes (input value minus one) unless 1 <= IER <= 4.
!

!       X,Y = Updated arrays of length N-1 containing nodal coordinates
!             (with elements K+1,...,N shifted up a position and thus
!             overwriting element K) unless 1 <= IER <= 4.  (N here denotes
!             the input value.)

!       LIST,LPTR,LEND,LNEW = Updated triangulation data structure
!                             reflecting the deletion unless IER .NE. 0.
!                             Note that the data structure may
!                             have been altered if IER >= 3.

!       LWK = Number of IWK columns required unless IER = 1 or IER = 3.

!       IWK = Indexes of the endpoints of the new arcs added unless LWK = 0
!             or 1 <= IER <= 4.  (Arcs are associated with columns, or
!             pairs of adjacent elements if IWK is declared as a
!             singly-subscripted array.)

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if K, NCC, N, or an LCC entry is outside its valid
!                     range or LWK < 0 on input.
!             IER = 2 if more space is required in IWK.  Refer to LWK.
!             IER = 3 if the triangulation data structure is invalid on input.
!             IER = 4 if K is an interior node with 4 or more neighbors, and
!                     the number of neighbors could not be reduced to 3 by
!                     swaps.  This could be caused by floating point errors
!                     with collinear nodes or by an invalid data structure.
!             IER = 5 if an error flag was returned by OPTIM.
!                     An error message is written to the standard output unit
!                     in this event.

!   Note that the deletion may result in all remaining nodes being collinear.
! This situation is not flagged.

! Modules required by DELNOD:  DELNB, LEFT, LSTPTR, NBCNT, OPTIM, SWAP, SWPTST

! Intrinsic function called by DELNOD:  ABS

!***********************************************************

INTEGER :: i, ierr, iwl, j, lccip1, lnw, lp, lp21, lpf,  &
           lph, lpl, lpl2, lpn, lwkl, n1, n2, nfrst, nit, nl, nn, nnb, nr
LOGICAL :: bdry
REAL    :: x1, x2, xl, xr, y1, y2, yl, yr

! Set N1 to K and NNB to the number of neighbors of N1 (plus
!   one if N1 is a boundary node), and test for errors.  LPF
!   and LPL are LIST indexes of the first and last neighbors
!   of N1, IWL is the number of IWK columns containing arcs,
!   and BDRY is TRUE iff N1 is a boundary node.

n1 = k
nn = n
IF (ncc < 0 .OR. n1 < 1 .OR. nn < 4  .OR. lwk < 0) GO TO 21
lccip1 = nn+1
DO  i = ncc,1,-1
  IF (lccip1-lcc(i) < 3) GO TO 21
  lccip1 = lcc(i)
END DO
IF (n1 >= lccip1) GO TO 21
lpl = lend(n1)
lpf = lptr(lpl)
nnb = nbcnt(lpl,lptr)
bdry = list(lpl) < 0
IF (bdry) nnb = nnb + 1
IF (nnb < 3) GO TO 23
lwkl = lwk
lwk = nnb - 3
IF (lwkl < lwk) GO TO 22
iwl = 0
IF (nnb == 3) GO TO 5

! Initialize for loop on arcs N1-N2 for neighbors N2 of N1, beginning with the
!   second neighbor.  NR and NL are the neighbors preceding and following N2,
!   respectively, and LP indexes NL.  The loop is exited when all possible
!   swaps have been applied to arcs incident on N1.  If N1 is interior,
!   the number of neighbors will be reduced to 3.

x1 = x(n1)
y1 = y(n1)
nfrst = list(lpf)
nr = nfrst
xr = x(nr)
yr = y(nr)
lp = lptr(lpf)
n2 = list(lp)
x2 = x(n2)
y2 = y(n2)
lp = lptr(lp)

! Top of loop:  set NL to the neighbor following N2.

2 nl = ABS(list(lp))
IF (nl == nfrst .AND. bdry) GO TO 5
xl = x(nl)
yl = y(nl)

!   Test for a convex quadrilateral.  To avoid an incorrect test caused by
!     collinearity, use the fact that if N1 is a boundary node, then
!     N1 LEFT NR->NL and if N2 is a boundary node, then N2 LEFT NL->NR.

lpl2 = lend(n2)
IF ( (bdry .OR. left(xr, yr, xl, yl, x1, y1)) .AND. (list(lpl2) < 0 .OR.  &
     left(xl, yl, xr, yr, x2, y2)) ) GO TO 3

!   Nonconvex quadrilateral -- no swap is possible.

nr = n2
xr = x2
yr = y2
GO TO 4

!   The quadrilateral defined by adjacent triangles (N1,N2,NL) and (N2,N1,NR)
!     is convex.  Swap in NL-NR and store it in IWK.  Indexes larger than N1
!     must be decremented since N1 will be deleted from X and Y.

3 CALL swap (nl, nr, n1, n2, list, lptr, lend, lp21)
iwl = iwl + 1
IF (nl <= n1) THEN
  iwk(1,iwl) = nl
ELSE
  iwk(1,iwl) = nl - 1
END IF
IF (nr <= n1) THEN
  iwk(2,iwl) = nr
ELSE
  iwk(2,iwl) = nr - 1
END IF

!   Recompute the LIST indexes LPL,LP and decrement NNB.

lpl = lend(n1)
nnb = nnb - 1
IF (nnb == 3) GO TO 5
lp = lstptr(lpl, nl, list, lptr)
IF (nr == nfrst) GO TO 4

!   NR is not the first neighbor of N1.
!     Back up and test N1-NR for a swap again:  Set N2 to NR and NR to the
!     previous neighbor of N1 -- the neighbor of NR which follows N1.
!     LP21 points to NL as a neighbor of NR.

n2 = nr
x2 = xr
y2 = yr
lp21 = lptr(lp21)
lp21 = lptr(lp21)
nr = ABS(list(lp21))
xr = x(nr)
yr = y(nr)
GO TO 2

!   Bottom of loop -- test for invalid termination.

4 IF (n2 == nfrst) GO TO 24
n2 = nl
x2 = xl
y2 = yl
lp = lptr(lp)
GO TO 2

! Delete N1 from the adjacency list of N2 for all neighbors N2 of N1.
!   LPL points to the last neighbor of N1.
!   LNEW is stored in local variable LNW.

5 lp = lpl
lnw = lnew

! Loop on neighbors N2 of N1, beginning with the first.

6 lp = lptr(lp)
n2 = ABS(list(lp))
CALL delnb (n2,n1,n, list,lptr,lend,lnw, lph)
IF (lph < 0) GO TO 23

!   LP and LPL may require alteration.

IF (lpl == lnw) lpl = lph
IF (lp == lnw) lp = lph
IF (lp /= lpl) GO TO 6

! Delete N1 from X, Y, and LEND, and remove its adjacency list from LIST and
!   LPTR.  LIST entries (nodal indexes) which are larger than N1 must be
!   decremented.

nn = nn - 1
IF (n1 > nn) GO TO 9
DO  i = n1,nn
  x(i) = x(i+1)
  y(i) = y(i+1)
  lend(i) = lend(i+1)
END DO

DO  i = 1,lnw-1
  IF (list(i) > n1) list(i) = list(i) - 1
  IF (list(i) < -n1) list(i) = list(i) + 1
END DO

!   For LPN = first to last neighbors of N1, delete the preceding neighbor
!     (indexed by LP).

!   Each empty LIST,LPTR location LP is filled in with the values at LNW-1,
!     and LNW is decremented.  All pointers (including those in LPTR and LEND)
!     with value LNW-1 must be changed to LP.

!  LPL points to the last neighbor of N1.

9 IF (bdry) nnb = nnb - 1
lpn = lpl
DO  j = 1,nnb
  lnw = lnw - 1
  lp = lpn
  lpn = lptr(lp)
  list(lp) = list(lnw)
  lptr(lp) = lptr(lnw)
  IF (lptr(lpn) == lnw) lptr(lpn) = lp
  IF (lpn == lnw) lpn = lp
  DO  i = nn,1,-1
    IF (lend(i) == lnw) THEN
      lend(i) = lp
      GO TO 11
    END IF
  END DO
  
  11   DO  i = lnw-1,1,-1
    IF (lptr(i) == lnw) lptr(i) = lp
  END DO
END DO

! Decrement LCC entries.

DO  i = 1,ncc
  lcc(i) = lcc(i) - 1
END DO

! Update N and LNEW, and optimize the patch of triangles
!   containing K (on input) by applying swaps to the arcs in IWK.

n = nn
lnew = lnw
IF (iwl > 0) THEN
  nit = 4*iwl
  CALL optim (x, y, iwl, list, lptr, lend, nit, iwk, ierr)
  IF (ierr /= 0) GO TO 25
END IF

! Successful termination.

ier = 0
RETURN

! Invalid input parameter.

21 ier = 1
RETURN

! Insufficient space reserved for IWK.

22 ier = 2
RETURN

! Invalid triangulation data structure.  NNB < 3 on input or
!   N2 is a neighbor of N1 but N1 is not a neighbor of N2.

23 ier = 3
RETURN

! K is an interior node with 4 or more neighbors, but the
!   number of neighbors could not be reduced.

24 ier = 4
RETURN

! Error flag returned by OPTIM.

25 ier = 5
WRITE (*,100) nit, ierr
RETURN
100 FORMAT (//t6, '*** Error in OPTIM:  NIT = ', i4, ', IER = ', i1, ' ***'/)
END SUBROUTINE delnod



SUBROUTINE edge (in1, in2, x, y, lwk, list, lptr, lend, ier)

! N.B. Argument IWK has been removed.

INTEGER, INTENT(IN)      :: in1
INTEGER, INTENT(IN)      :: in2
REAL, INTENT(IN)         :: x(:)
REAL, INTENT(IN)         :: y(:)
INTEGER, INTENT(IN OUT)  :: lwk
INTEGER, INTENT(IN OUT)  :: list(:)
INTEGER, INTENT(IN OUT)  :: lptr(:)
INTEGER, INTENT(IN OUT)  :: lend(:)
INTEGER, INTENT(OUT)     :: ier

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   06/23/98

!   Given a triangulation of N nodes and a pair of nodal indexes IN1 and IN2,
! this routine swaps arcs as necessary to force IN1 and IN2 to be adjacent.
! Only arcs which intersect IN1-IN2 are swapped out.  If a Delaunay
! triangulation is input, the resulting triangulation is as close as possible
! to a Delaunay triangulation in the sense that all arcs other than IN1-IN2
! are locally optimal.

!   A sequence of calls to EDGE may be used to force the presence of a set of
! edges defining the boundary of a non-convex and/or multiply connected region
! (refer to Subroutine ADDCST), or to introduce barriers into the
! triangulation.  Note that Subroutine GETNP will not necessarily return
! closest nodes if the triangulation has been constrained by a call to EDGE.
! However, this is appropriate in some applications, such as triangle-based
! interpolation on a nonconvex domain.

! On input:

!       IN1,IN2 = Indexes (of X and Y) in the range 1 to N
!                 defining a pair of nodes to be connected by an arc.

!       X,Y = Arrays of length N containing the Cartesian coordinates of the
!             nodes.

! The above parameters are not altered by this routine.

!       LWK = Number of columns reserved for IWK.  This must be at least NI --
!             the number of arcs which intersect IN1-IN2.
!             (NI is bounded by N-3.)

!       IWK = Integer work array of length at least 2*LWK.

!       LIST,LPTR,LEND = Data structure defining the triangulation.
!                        Refer to Subroutine TRMESH.

! On output:

!       LWK = Number of arcs which intersect IN1-IN2 (but not more than the
!             input value of LWK) unless IER = 1 or IER = 3.  LWK = 0 if and
!             only if IN1 and IN2 were adjacent (or LWK=0) on input.

!       IWK = Array containing the indexes of the endpoints of the new arcs
!             other than IN1-IN2 unless IER > 0 or LWK = 0.  New arcs to the
!             left of IN2-IN1 are stored in the first K-1 columns
!             (left portion of IWK), column K contains zeros, and new arcs to
!             the right of IN2-IN1 occupy columns K+1,...,LWK.  (K can be
!             determined by searching IWK for the zeros.)

!       LIST,LPTR,LEND = Data structure updated if necessary to reflect the
!                        presence of an arc connecting IN1 and IN2 unless
!                        IER.NE. 0.  The data structure has been altered if
!                        IER = 4.

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if IN1 < 1, IN2 < 1, IN1 = IN2, or LWK < 0 on input.
!             IER = 2 if more space is required in IWK.
!             IER = 3 if IN1 and IN2 could not be connected due to
!                     either an invalid data structure or
!                     collinear nodes (and floating point error).
!             IER = 4 if an error flag was returned by OPTIM.

!   An error message is written to the standard output unit
! in the case of IER = 3 or IER = 4.

! Modules required by EDGE:  LEFT, LSTPTR, OPTIM, SWAP, SWPTST

! Intrinsic function called by EDGE:  ABS

!***********************************************************

INTEGER :: i, ierr, iwc, iwcp1, iwend, iwf, iwk(2,lwk), iwl, lft, lp,  &
           lpl, lp21, next, nit, nl, nr, n0, n1, n2, n1frst, n1lst
REAL    :: dx, dy, x0, y0, x1, y1, x2, y2

! Local parameters:

! DX,DY =   Components of arc N1-N2
! I =       DO-loop index and column index for IWK
! IERR =    Error flag returned by Subroutine OPTIM
! IWC =     IWK index between IWF and IWL -- NL->NR is
!             stored in IWK(1,IWC)->IWK(2,IWC)
! IWCP1 =   IWC + 1
! IWEND =   Input or output value of LWK
! IWF =     IWK (column) index of the first (leftmost) arc
!             which intersects IN1->IN2
! IWL =     IWK (column) index of the last (rightmost) are
!             which intersects IN1->IN2
! LFT =     Flag used to determine if a swap results in the new arc
!             intersecting IN1-IN2 -- LFT = 0 iff N0 = IN1, LFT = -1 implies
!             N0 LEFT IN1->IN2, and LFT = 1 implies N0 LEFT IN2->IN1
! LP21 =    Unused parameter returned by SWAP
! LP =      List pointer (index) for LIST and LPTR
! LPL =     Pointer to the last neighbor of IN1 or NL
! N0 =      Neighbor of N1 or node opposite NR->NL
! N1,N2 =   Local copies of IN1 and IN2
! N1FRST =  First neighbor of IN1
! N1LST =   (Signed) last neighbor of IN1
! NEXT =    Node opposite NL->NR
! NIT =     Flag or number of iterations employed by OPTIM
! NL,NR =   Endpoints of an arc which intersects IN1-IN2 with NL LEFT IN1->IN2
! X0,Y0 =   Coordinates of N0
! X1,Y1 =   Coordinates of IN1
! X2,Y2 =   Coordinates of IN2


! Store IN1, IN2, and LWK in local variables and test for errors.

n1 = in1
n2 = in2
iwend = lwk
IF (n1 < 1 .OR. n2 < 1 .OR. n1 == n2  .OR. iwend < 0) GO TO 31

! Test for N2 as a neighbor of N1.  LPL points to the last neighbor of N1.

lpl = lend(n1)
n0 = ABS(list(lpl))
lp = lpl
1 IF (n0 == n2) GO TO 30
lp = lptr(lp)
n0 = list(lp)
IF (lp /= lpl) GO TO 1

! Initialize parameters.

iwl = 0
nit = 0

! Store the coordinates of N1 and N2.

2 x1 = x(n1)
y1 = y(n1)
x2 = x(n2)
y2 = y(n2)

! Set NR and NL to adjacent neighbors of N1 such that NR LEFT N2->N1 and
!   NL LEFT N1->N2, (NR Forward N1->N2 or NL Forward N1->N2), and
!   (NR Forward N2->N1 or NL Forward N2->N1).

!   Initialization:  Set N1FRST and N1LST to the first and (signed) last
!     neighbors of N1, respectively, and initialize NL to N1FRST.

lpl = lend(n1)
n1lst = list(lpl)
lp = lptr(lpl)
n1frst = list(lp)
nl = n1frst
IF (n1lst < 0) GO TO 4

!   N1 is an interior node.  Set NL to the first candidate
!     for NR (NL LEFT N2->N1).

3 IF ( left(x2,y2,x1,y1,x(nl),y(nl)) ) GO TO 4
lp = lptr(lp)
nl = list(lp)
IF (nl /= n1frst) GO TO 3

!   All neighbors of N1 are strictly left of N1->N2.

GO TO 5

!   NL = LIST(LP) LEFT N2->N1.  Set NR to NL and NL to the
!     following neighbor of N1.

4 nr = nl
lp = lptr(lp)
nl = ABS(list(lp))
IF ( left(x1,y1,x2,y2,x(nl),y(nl)) ) THEN
  
!   NL LEFT N1->N2 and NR LEFT N2->N1.  The Forward tests are employed to
!     avoid an error associated with collinear nodes.
  
  dx = x2-x1
  dy = y2-y1
  IF ((dx*(x(nl)-x1)+dy*(y(nl)-y1) >= 0. .OR. &
      dx*(x(nr)-x1)+dy*(y(nr)-y1) >= 0.) .AND. &
      (dx*(x(nl)-x2)+dy*(y(nl)-y2) <= 0. .OR. &
      dx*(x(nr)-x2)+dy*(y(nr)-y2) <= 0.)) GO TO 6
  
!   NL-NR does not intersect N1-N2.  However, there is another
!     candidate for the first arc if NL lies on the line N1-N2.
  
  IF ( .NOT. left(x2,y2,x1,y1,x(nl),y(nl)) ) GO TO 5
END IF

!   Bottom of loop.

IF (nl /= n1frst) GO TO 4

! Either the triangulation is invalid or N1-N2 lies on the convex hull
!   boundary and an edge NR->NL (opposite N1 and intersecting N1-N2)
!   was not found due to floating point error.
!   Try interchanging N1 and N2 -- NIT > 0 iff this has already been done.

5 IF (nit > 0) GO TO 33
nit = 1
n1 = n2
n2 = in1
GO TO 2

! Store the ordered sequence of intersecting edges NL->NR in
!   IWK(1,IWL)->IWK(2,IWL).

6 iwl = iwl + 1
IF (iwl > iwend) GO TO 32
iwk(1,iwl) = nl
iwk(2,iwl) = nr

!   Set NEXT to the neighbor of NL which follows NR.

lpl = lend(nl)
lp = lptr(lpl)

!   Find NR as a neighbor of NL.  The search begins with the first neighbor.

7 IF (list(lp) == nr) GO TO 8
lp = lptr(lp)
IF (lp /= lpl) GO TO 7

!   NR must be the last neighbor, and NL->NR cannot be a boundary edge.

IF (list(lp) /= nr) GO TO 33

!   Set NEXT to the neighbor following NR, and test for
!     termination of the store loop.

8 lp = lptr(lp)
next = ABS(list(lp))
IF (next == n2) GO TO 9

!   Set NL or NR to NEXT.

IF ( left(x1,y1,x2,y2,x(next),y(next)) ) THEN
  nl = next
ELSE
  nr = next
END IF
GO TO 6

! IWL is the number of arcs which intersect N1-N2.
!   Store LWK.

9 lwk = iwl
iwend = iwl

! Initialize for edge swapping loop -- all possible swaps are applied
!   (even if the new arc again intersects N1-N2), arcs to the left of N1->N2
!   are stored in the left portion of IWK, and arcs to the right are stored
!   in the right portion.  IWF and IWL index the first and last
!   intersecting arcs.

iwf = 1

! Top of loop -- set N0 to N1 and NL->NR to the first edge.
!   IWC points to the arc currently being processed.
!   LFT <= 0 iff N0 LEFT N1->N2.

10 lft = 0
n0 = n1
x0 = x1
y0 = y1
nl = iwk(1,iwf)
nr = iwk(2,iwf)
iwc = iwf

!   Set NEXT to the node opposite NL->NR unless IWC is the last arc.

11 IF (iwc == iwl) GO TO 21
iwcp1 = iwc + 1
next = iwk(1,iwcp1)
IF (next /= nl) GO TO 16
next = iwk(2,iwcp1)

!   NEXT RIGHT N1->N2 and IWC < IWL.  Test for a possible swap.

IF ( .NOT. left(x0,y0,x(nr),y(nr),x(next),y(next)) ) GO TO 14
IF (lft >= 0) GO TO 12
IF ( .NOT. left(x(nl),y(nl),x0,y0,x(next),y(next)) ) GO TO 14

!   Replace NL->NR with N0->NEXT.

CALL swap (next,n0,nl,nr, list,lptr,lend, lp21)
iwk(1,iwc) = n0
iwk(2,iwc) = next
GO TO 15

!   Swap NL-NR for N0-NEXT, shift columns IWC+1,...,IWL to
!     the left, and store N0-NEXT in the right portion of IWK.

12 CALL swap (next,n0,nl,nr, list,lptr,lend, lp21)
DO  i = iwcp1,iwl
  iwk(1,i-1) = iwk(1,i)
  iwk(2,i-1) = iwk(2,i)
END DO
iwk(1,iwl) = n0
iwk(2,iwl) = next
iwl = iwl - 1
nr = next
GO TO 11

!   A swap is not possible.  Set N0 to NR.

14 n0 = nr
x0 = x(n0)
y0 = y(n0)
lft = 1

!   Advance to the next arc.

15 nr = next
iwc = iwc + 1
GO TO 11

!   NEXT LEFT N1->N2, NEXT .NE. N2, and IWC < IWL.
!     Test for a possible swap.

16 IF ( .NOT. left(x(nl),y(nl),x0,y0,x(next),y(next)) ) GO TO 19
IF (lft <= 0) GO TO 17
IF ( .NOT. left(x0,y0,x(nr),y(nr),x(next),y(next)) ) GO TO 19

!   Replace NL->NR with NEXT->N0.

CALL swap (next, n0, nl, nr, list, lptr, lend, lp21)
iwk(1,iwc) = next
iwk(2,iwc) = n0
GO TO 20

!   Swap NL-NR for N0-NEXT, shift columns IWF,...,IWC-1 to
!     the right, and store N0-NEXT in the left portion of
!     IWK.

17 CALL swap (next, n0, nl, nr, list, lptr, lend, lp21)
DO  i = iwc-1,iwf,-1
  iwk(1,i+1) = iwk(1,i)
  iwk(2,i+1) = iwk(2,i)
END DO
iwk(1,iwf) = n0
iwk(2,iwf) = next
iwf = iwf + 1
GO TO 20

!   A swap is not possible.  Set N0 to NL.

19 n0 = nl
x0 = x(n0)
y0 = y(n0)
lft = -1

!   Advance to the next arc.

20 nl = next
iwc = iwc + 1
GO TO 11

!   N2 is opposite NL->NR (IWC = IWL).

21 IF (n0 == n1) GO TO 24
IF (lft < 0) GO TO 22

!   N0 RIGHT N1->N2.  Test for a possible swap.

IF ( .NOT. left(x0, y0, x(nr), y(nr), x2, y2) ) GO TO 10

!   Swap NL-NR for N0-N2 and store N0-N2 in the right portion of IWK.

CALL swap (n2, n0, nl, nr, list, lptr, lend, lp21)
iwk(1,iwl) = n0
iwk(2,iwl) = n2
iwl = iwl - 1
GO TO 10

!   N0 LEFT N1->N2.  Test for a possible swap.

22 IF ( .NOT. left(x(nl), y(nl), x0, y0, x2, y2) ) GO TO 10

!   Swap NL-NR for N0-N2, shift columns IWF,...,IWL-1 to the
!     right, and store N0-N2 in the left portion of IWK.

CALL swap (n2, n0, nl, nr, list, lptr, lend, lp21)
i = iwl
23 iwk(1,i) = iwk(1,i-1)
iwk(2,i) = iwk(2,i-1)
i = i - 1
IF (i > iwf) GO TO 23
iwk(1,iwf) = n0
iwk(2,iwf) = n2
iwf = iwf + 1
GO TO 10

! IWF = IWC = IWL.  Swap out the last arc for N1-N2 and store zeros in IWK.

24 CALL swap (n2, n1, nl, nr, list, lptr, lend, lp21)
iwk(1,iwc) = 0
iwk(2,iwc) = 0

! Optimization procedure --

IF (iwc > 1) THEN
  
!   Optimize the set of new arcs to the left of IN1->IN2.
  
  nit = 3*(iwc-1)
  CALL optim (x, y, iwc-1, list, lptr, lend, nit, iwk, ierr)
  IF (ierr /= 0) GO TO 34
END IF
IF (iwc < iwend) THEN
  
!   Optimize the set of new arcs to the right of IN1->IN2.
  
  nit = 3*(iwend-iwc)
  CALL optim (x, y, iwend-iwc, list, lptr, lend, nit, iwk(1:,iwc+1:), ierr)
  IF (ierr /= 0) GO TO 34
END IF

! Successful termination.

ier = 0
RETURN

! IN1 and IN2 were adjacent on input.

30 ier = 0
RETURN

! Invalid input parameter.

31 ier = 1
RETURN

! Insufficient space reserved for IWK.

32 ier = 2
RETURN

! Invalid triangulation data structure or collinear nodes
!   on convex hull boundary.

33 ier = 3
WRITE (*,130) in1, in2
130 FORMAT (//t6, '*** Error in EDGE:  Invalid triangulation',  &
            ' or null triangles on boundary'/ t10, 'IN1 =', i4, ', IN2 =', i4/)
RETURN

! Error flag returned by OPTIM.

34 ier = 4
WRITE (*,140) nit, ierr
140 FORMAT (//t6, '*** Error in OPTIM:  NIT = ', i4, ', IER = ', i1, ' ***'/)
RETURN
END SUBROUTINE edge



SUBROUTINE getnp (ncc, lcc, n, x, y, list, lptr, lend, l, npts, ds, ier)

INTEGER, INTENT(IN)      :: ncc
INTEGER, INTENT(IN)      :: lcc(:)
INTEGER, INTENT(IN)      :: n
REAL, INTENT(IN)         :: x(:)
REAL, INTENT(IN)         :: y(:)
INTEGER, INTENT(IN OUT)  :: list(:)
INTEGER, INTENT(IN)      :: lptr(:)
INTEGER, INTENT(OUT)     :: lend(:)
INTEGER, INTENT(IN)      :: l
INTEGER, INTENT(IN OUT)  :: npts(:)
REAL, INTENT(IN OUT)     :: ds(:)
INTEGER, INTENT(OUT)     :: ier

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   11/12/94

!   Given a triangulation of N nodes and an array NPTS con-
! taining the indexes of L-1 nodes ordered by distance from
! NPTS(1), this subroutine sets NPTS(L) to the index of the
! next node in the sequence -- the node, other than NPTS(1),
! ...,NPTS(L-1), which is closest to NPTS(1).  Thus, the
! ordered sequence of K closest nodes to N1 (including N1)
! may be determined by K-1 calls to GETNP with NPTS(1) = N1
! and L = 2,3,...,K for K >= 2.  Note that NPTS must in-
! clude constraint nodes as well as non-constraint nodes.
! Thus, a sequence of K1 closest non-constraint nodes to N1
! must be obtained as a subset of the closest K2 nodes to N1
! for some K2 >= K1.

!   The terms closest and distance have special definitions
! when constraint nodes are present in the triangulation.
! Nodes N1 and N2 are said to be visible from each other if
! and only if the line segment N1-N2 intersects no con-
! straint arc (except possibly itself) and is not an interi-
! or constraint arc (arc whose interior lies in a constraint
! region).  A path from N1 to N2 is an ordered sequence of
! nodes, with N1 first and N2 last, such that adjacent path
! elements are visible from each other.  The path length is
! the sum of the Euclidean distances between adjacent path
! nodes.  Finally, the distance from N1 to N2 is defined to
! be the length of the shortest path from N1 to N2.

!   The algorithm uses the property of a Delaunay triangula-
! tion that the K-th closest node to N1 is a neighbor of one
! of the K-1 closest nodes to N1.  With the definition of
! distance used here, this property holds when constraints
! are present as long as non-constraint arcs are locally optimal.


! On input:

!       NCC = Number of constraints.  NCC >= 0.

!       LCC = List of constraint curve starting indexes (or dummy array of
!             length 1 if NCC = 0).  Refer to Subroutine ADDCST.

!       N = Number of nodes in the triangulation.  N >= 3.

!       X,Y = Arrays of length N containing the coordinates
!             of the nodes with non-constraint nodes in the
!             first LCC(1)-1 locations if NCC > 0.

!       LIST,LPTR,LEND = Triangulation data structure.  Refer to
!                        Subroutine TRMESH.

!       L = Number of nodes in the sequence on output.  2 <= L <= N.

!       NPTS = Array of length >= L containing the indexes of the L-1 closest
!              nodes to NPTS(1) in the first L-1 locations.

!       DS = Array of length >= L containing the distance
!            (defined above) between NPTS(1) and NPTS(I) in
!            the I-th position for I = 1,...,L-1.  Thus, DS(1) = 0.

! Input parameters other than NPTS(L) and DS(L) are not
!   altered by this routine.

! On output:

!       NPTS = Array updated with the index of the L-th closest node
!              to NPTS(1) in position L unless IER .NE. 0.

!       DS = Array updated with the distance between NPTS(1)
!            and NPTS(L) in position L unless IER .NE. 0.

!       IER = Error indicator:
!             IER =  0 if no errors were encountered.
!             IER = -1 if NCC, N, L, or an LCC entry is
!                      outside its valid range on input.
!             IER =  K if NPTS(K) is not a valid index in the range 1 to N.

! Module required by GETNP:  INTSEC

! Intrinsic functions called by GETNP:  ABS, MIN, SQRT

!***********************************************************

INTEGER :: i, ifrst, ilast, j, k, km1, lcc1, lm1, lp,  &
           lpcl, lpk, lpkl, n1, nc, nf1, nf2, nj, nk, nkbak, nkfor, nl, nn
LOGICAL :: isw, vis, ncf, njf, skip, sksav, lft1, lft2, lft12
REAL    :: dc, dl, x1, xc, xj, xk, y1, yc, yj, yk

! Store parameters in local variables and test for errors.
!   LCC1 indexes the first constraint node.

ier = -1
nn = n
lcc1 = nn+1
lm1 = l-1
IF (ncc < 0 .OR. lm1 < 1 .OR. lm1 >= nn) RETURN
IF (ncc == 0) THEN
  IF (nn < 3) RETURN
ELSE
  DO  i = ncc,1,-1
    IF (lcc1 - lcc(i) < 3) RETURN
    lcc1 = lcc(i)
  END DO
  IF (lcc1 < 1) RETURN
END IF

! Test for an invalid index in NPTS.

DO  k = 1,lm1
  nk = npts(k)
  IF (nk < 1 .OR. nk > nn) THEN
    ier = k
    RETURN
  END IF
END DO

! Store N1 = NPTS(1) and mark the elements of NPTS.

n1 = npts(1)
x1 = x(n1)
y1 = y(n1)
DO  k = 1,lm1
  nk = npts(k)
  lend(nk) = -lend(nk)
END DO

! Candidates NC for NL = NPTS(L) are the unmarked visible
!   neighbors of nodes NK in NPTS.  ISW is an initialization
!   switch set to .TRUE. when NL and its distance DL from N1
!   have been initialized with the first candidate encountered.

isw = .false.
dl = 0.

! Loop on marked nodes NK = NPTS(K).  LPKL indexes the last
!   neighbor of NK in LIST.

DO  k = 1,lm1
  km1 = k - 1
  nk = npts(k)
  xk = x(nk)
  yk = y(nk)
  lpkl = -lend(nk)
  nkfor = 0
  nkbak = 0
  vis = .true.
  IF (nk >= lcc1) THEN
    
!   NK is a constraint node.  Set NKFOR and NKBAK to the constraint nodes
!     which follow and precede NK.  IFRST and ILAST are set to the first and
!     last nodes in the constraint containing NK.
    
    ifrst = nn + 1
    DO  i = ncc,1,-1
      ilast = ifrst - 1
      ifrst = lcc(i)
      IF (nk >= ifrst) GO TO 5
    END DO
    
    5 IF (nk < ilast) THEN
      nkfor = nk + 1
    ELSE
      nkfor = ifrst
    END IF
    IF (nk > ifrst) THEN
      nkbak = nk - 1
    ELSE
      nkbak = ilast
    END IF
    
!   Initialize VIS to TRUE iff NKFOR precedes NKBAK in the adjacency list for
!     NK -- the first neighbor is visible and is not NKBAK.
    
    lpk = lpkl
    6 lpk = lptr(lpk)
    nc = ABS(list(lpk))
    IF (nc /= nkfor .AND. nc /= nkbak) GO TO 6
    vis = nc == nkfor
  END IF
  
! Loop on neighbors NC of NK, bypassing marked and nonvisible neighbors.
  
  lpk = lpkl
  7 lpk = lptr(lpk)
  nc = ABS(list(lpk))
  IF (nc == nkbak) vis = .true.
  
!   VIS = .FALSE. iff NK-NC is an interior constraint arc
!     (NK is a constraint node and NC lies strictly between NKFOR and NKBAK).
  
  IF (.NOT. vis) GO TO 15
  IF (nc == nkfor) vis = .false.
  IF (lend(nc) < 0) GO TO 15
  
! Initialize distance DC between N1 and NC to Euclidean distance.
  
  xc = x(nc)
  yc = y(nc)
  dc = SQRT((xc-x1)*(xc-x1) + (yc-y1)*(yc-y1))
  IF (isw .AND. dc >= dl) GO TO 15
  IF (k == 1) GO TO 14
  
! K >= 2.  Store the pointer LPCL to the last neighbor of NC.
  
  lpcl = lend(nc)
  
! Set DC to the length of the shortest path from N1 to NC
!   which has not previously been encountered and which is
!   a viable candidate for the shortest path from N1 to NL.
!   This is Euclidean distance iff NC is visible from N1.
!   Since the shortest path from N1 to NL contains only ele-
!   ments of NPTS which are constraint nodes (in addition to
!   N1 and NL), only these need be considered for the path
!   from N1 to NC.  Thus, for distance function D(A,B) and
!   J = 1,...,K, DC = min(D(N1,NJ) + D(NJ,NC)) over con-
!   straint nodes NJ = NPTS(J) which are visible from NC.
  
  DO  j = 1,km1
    nj = npts(j)
    IF (j > 1 .AND. nj < lcc1) CYCLE
    
! If NC is a visible neighbor of NJ, a path from N1 to NC
!   containing NJ has already been considered.  Thus, NJ may
!   be bypassed if it is adjacent to NC.
    
    lp = lpcl
    8 lp = lptr(lp)
    IF ( nj == ABS(list(lp)) ) GO TO 12
    IF (lp /= lpcl) GO TO 8
    
! NJ is a constraint node (unless J=1) not adjacent to NC,
!   and is visible from NC iff NJ-NC is not intersected by
!   a constraint arc.  Loop on constraints I in reverse order --
    
    xj = x(nj)
    yj = y(nj)
    ifrst = nn+1
    DO  i = ncc,1,-1
      ilast = ifrst - 1
      ifrst = lcc(i)
      nf1 = ilast
      ncf = nf1 == nc
      njf = nf1 == nj
      skip = ncf .OR. njf
      
! Loop on boundary constraint arcs NF1-NF2 which contain neither NC nor NJ.
!   NCF and NJF are TRUE iff NC (or NJ) has been encountered in the
!   constraint, and SKIP = .TRUE. iff NF1 = NC or NF1 = NJ.
      
      DO  nf2 = ifrst,ilast
        IF (nf2 == nc) ncf = .true.
        IF (nf2 == nj) njf = .true.
        sksav = skip
        skip = nf2 == nc .OR. nf2 == nj
        
!   The last constraint arc in the constraint need not be
!     tested if none of the arcs have been skipped.
        
        IF ( sksav .OR. skip  .OR. (nf2 == ilast .AND. &
            .NOT. ncf .AND. .NOT. njf) ) GO TO 9
        IF ( intsec(x(nf1),y(nf1),x(nf2),y(nf2), xc,yc,xj,yj) ) GO TO 12
        9 nf1 = nf2
      END DO
      IF (.NOT. ncf .OR. .NOT. njf) CYCLE
      
! NC and NJ are constraint nodes in the same constraint.
!   NC-NJ is intersected by an interior constraint arc iff
!   1)  NC LEFT NF2->NF1 and (NJ LEFT NF1->NC and NJ LEFT NC->NF2) or
!   2)  NC .NOT. LEFT NF2->NF1 and (NJ LEFT NF1->NC or NJ LEFT NC->NF2),
!   where NF1, NC, NF2 are consecutive constraint nodes.
      
      IF (nc /= ifrst) THEN
        nf1 = nc - 1
      ELSE
        nf1 = ilast
      END IF
      IF (nc /= ilast) THEN
        nf2 = nc + 1
      ELSE
        nf2 = ifrst
      END IF
      lft1 = (xc-x(nf1))*(yj-y(nf1)) >= (xj-x(nf1))*(yc-y(nf1))
      lft2 = (x(nf2)-xc)*(yj-yc) >= (xj-xc)*(y(nf2)-yc)
      lft12 = (x(nf1)-x(nf2))*(yc-y(nf2)) >= (xc-x(nf2))*(y(nf1)-y(nf2))
      IF ( (lft1 .AND. lft2) .OR. (.NOT. lft12  &
          .AND. (lft1 .OR. lft2)) ) GO TO 12
    END DO
    
! NJ is visible from NC.  Exit the loop with DC = Euclidean distance if J = 1.
    
    IF (j == 1) GO TO 14
    dc = MIN(dc,ds(j) + SQRT((xc-xj)*(xc-xj) + (yc-yj)*(yc-yj)))
    CYCLE
    
! NJ is not visible from NC or is adjacent to NC.  Initial-
!   ize DC with D(N1,NK) + D(NK,NC) if J = 1.
    
    12 IF (j == 1) dc = ds(k) + SQRT((xc-xk)*(xc-xk) + (yc-yk)*(yc-yk))
  END DO
  
! Compare DC with DL.
  
  IF (isw .AND. dc >= dl) GO TO 15
  
! The first (or a closer) candidate for NL has been encountered.
  
  14 nl = nc
  dl = dc
  isw = .true.
  15 IF (lpk /= lpkl) GO TO 7
END DO

! Unmark the elements of NPTS and store NL and DL.

DO  k = 1,lm1
  nk = npts(k)
  lend(nk) = -lend(nk)
END DO
npts(l) = nl
ds(l) = dl
ier = 0
RETURN
END SUBROUTINE getnp



FUNCTION indxcc (ncc, lcc, n, list, lend) RESULT(fn_val)

INTEGER, INTENT(IN)  :: ncc
INTEGER, INTENT(IN)  :: lcc(:)
INTEGER, INTENT(IN)  :: n
INTEGER, INTENT(IN)  :: list(:)
INTEGER, INTENT(IN)  :: lend(:)
INTEGER              :: fn_val

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   08/25/91

!   Given a constrained Delaunay triangulation, this func-
! tion returns the index, if any, of an exterior constraint
! curve (an unbounded constraint region).  An exterior con-
! straint curve is assumed to be present if and only if the
! clockwise-ordered sequence of boundary nodes is a subse-
! quence of a constraint node sequence.  The triangulation
! adjacencies corresponding to constraint edges may or may
! not have been forced by a call to ADDCST, and the con-
! straint region may or may not be valid (contain no nodes).


! On input:

!       NCC = Number of constraints.  NCC >= 0.

!       LCC = List of constraint curve starting indexes (or dummy array of
!             length 1 if NCC = 0).  Refer to Subroutine ADDCST.

!       N = Number of nodes in the triangulation.  N >= 3.

!       LIST,LEND = Data structure defining the triangulation.
!                   Refer to Subroutine TRMESH.

!   Input parameters are not altered by this function.  Note
! that the parameters are not tested for validity.

! On output:

!       INDXCC = Index of the exterior constraint curve, if
!                present, or 0 otherwise.

! Modules required by INDXCC:  None

!***********************************************************

INTEGER :: i, ifrst, ilast, lp, n0, nst, nxt

fn_val = 0
IF (ncc < 1) RETURN

! Set N0 to the boundary node with smallest index.

n0 = 0
1 n0 = n0 + 1
lp = lend(n0)
IF (list(lp) > 0) GO TO 1

! Search in reverse order for the constraint I, if any, that contains N0.
!   IFRST and ILAST index the first and last nodes in constraint I.

i = ncc
ilast = n
2 ifrst = lcc(i)
IF (n0 >= ifrst) GO TO 3
IF (i == 1) RETURN
i = i - 1
ilast = ifrst - 1
GO TO 2

! N0 is in constraint I which indexes an exterior constraint curve iff the
!   clockwise-ordered sequence of boundary node indexes beginning with N0 is
!   increasing and bounded above by ILAST.

3 nst = n0

4 nxt = -list(lp)
IF (nxt == nst) GO TO 5
IF (nxt <= n0 .OR. nxt > ilast) RETURN
n0 = nxt
lp = lend(n0)
GO TO 4

! Constraint I contains the boundary node sequence as a subset.

5 fn_val = i
RETURN
END FUNCTION indxcc



SUBROUTINE insert (k, lp, list, lptr, lnew )

INTEGER, INTENT(IN)      :: k
INTEGER, INTENT(IN OUT)  :: lp
INTEGER, INTENT(OUT)     :: list(:)
INTEGER, INTENT(IN OUT)  :: lptr(:)
INTEGER, INTENT(IN OUT)  :: lnew

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88

!   This subroutine inserts K as a neighbor of N1 following
! N2, where LP is the LIST pointer of N2 as a neighbor of
! N1.  Note that, if N2 is the last neighbor of N1, K will
! become the first neighbor (even if N1 is a boundary node).


! On input:

!       K = Index of the node to be inserted.

!       LP = LIST pointer of N2 as a neighbor of N1.

! The above parameters are not altered by this routine.

!       LIST,LPTR,LNEW = Data structure defining the triangulation.
!                        Refer to Subroutine TRMESH.

! On output:

!       LIST,LPTR,LNEW = Data structure updated with the
!                        addition of node K.

! Modules required by INSERT:  None

!***********************************************************

INTEGER :: lsav

lsav = lptr(lp)
lptr(lp) = lnew
list(lnew) = k
lptr(lnew) = lsav
lnew = lnew + 1
RETURN
END SUBROUTINE insert



SUBROUTINE intadd (kk, i1, i2, i3, list, lptr, lend, lnew )

INTEGER, INTENT(IN)      :: kk
INTEGER, INTENT(IN)      :: i1
INTEGER, INTENT(IN)      :: i2
INTEGER, INTENT(IN)      :: i3
INTEGER, INTENT(IN OUT)  :: list(:)
INTEGER, INTENT(OUT)     :: lptr(:)
INTEGER, INTENT(OUT)     :: lend(:)
INTEGER, INTENT(IN OUT)  :: lnew

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/22/91

!   This subroutine adds an interior node to a triangulation
! of a set of points in the plane.  The data structure is
! updated with the insertion of node KK into the triangle
! whose vertices are I1, I2, and I3.  No optimization of the
! triangulation is performed.


! On input:

!       KK = Index of the node to be inserted.  KK >= 1
!            and KK must not be equal to I1, I2, or I3.

!       I1,I2,I3 = Indexes of the counterclockwise-ordered sequence of
!                  vertices of a triangle which contains node KK.

! The above parameters are not altered by this routine.

!       LIST,LPTR,LEND,LNEW = Data structure defining the triangulation.
!                             Refer to Subroutine TRMESH.  Triangle (I1,I2,I3)
!                             must be included in the triangulation.

! On output:

!       LIST,LPTR,LEND,LNEW = Data structure updated with
!                             the addition of node KK.  KK
!                             will be connected to nodes I1,
!                             I2, and I3.

! Modules required by INTADD:  INSERT, LSTPTR

!***********************************************************

INTEGER :: k, lp, n1, n2, n3

k = kk

! Initialization.

n1 = i1
n2 = i2
n3 = i3

! Add K as a neighbor of I1, I2, and I3.

lp = lstptr(lend(n1),n2,list,lptr)
CALL insert (k,lp,list,lptr,lnew)
lp = lstptr(lend(n2),n3,list,lptr)
CALL insert (k,lp,list,lptr,lnew)
lp = lstptr(lend(n3),n1,list,lptr)
CALL insert (k,lp,list,lptr,lnew)

! Add I1, I2, and I3 as neighbors of K.

list(lnew) = n1
list(lnew+1) = n2
list(lnew+2) = n3
lptr(lnew) = lnew + 1
lptr(lnew+1) = lnew + 2
lptr(lnew+2) = lnew
lend(k) = lnew + 2
lnew = lnew + 3
RETURN
END SUBROUTINE intadd



FUNCTION intsec (x1, y1, x2, y2, x3, y3, x4, y4) RESULT(fn_val)

REAL, INTENT(IN)  :: x1
REAL, INTENT(IN)  :: y1
REAL, INTENT(IN)  :: x2
REAL, INTENT(IN)  :: y2
REAL, INTENT(IN)  :: x3
REAL, INTENT(IN)  :: y3
REAL, INTENT(IN)  :: x4
REAL, INTENT(IN)  :: y4
LOGICAL           :: fn_val

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88

!   Given a pair of line segments P1-P2 and P3-P4, this
! function returns the value .TRUE. if and only if P1-P2
! shares one or more points with P3-P4.  The line segments
! include their endpoints, and the four points need not be
! distinct.  Thus, either line segment may consist of a
! single point, and the segments may meet in a V (which is
! treated as an intersection).  Note that an incorrect
! decision may result from floating point error if the four
! endpoints are nearly collinear.


! On input:

!       X1,Y1 = Coordinates of P1.

!       X2,Y2 = Coordinates of P2.

!       X3,Y3 = Coordinates of P3.

!       X4,Y4 = Coordinates of P4.

! Input parameters are not altered by this function.

! On output:

!       INTSEC = Logical value defined above.

! Modules required by INTSEC:  None

!***********************************************************

REAL :: a, b, d, dx12, dx31, dx34, dy12, dy31, dy34

! Test for overlap between the smallest rectangles that
!   contain the line segments and have sides parallel to the axes.

IF ((x1 < x3 .AND. x1 < x4 .AND. x2 < x3 .AND.  x2 < x4) .OR. &
      (x1 > x3 .AND. x1 > x4 .AND. x2 > x3 .AND.  x2 > x4) .OR. &
      (y1 < y3 .AND. y1 < y4 .AND. y2 < y3 .AND.  y2 < y4) .OR. &
      (y1 > y3 .AND. y1 > y4 .AND. y2 > y3 .AND.  y2 > y4)) THEN
  fn_val = .false.
  RETURN
END IF

! Compute A = P4-P3 X P1-P3, B = P2-P1 X P1-P3, and
!   D = P2-P1 X P4-P3 (Z components).

dx12 = x2 - x1
dy12 = y2 - y1
dx34 = x4 - x3
dy34 = y4 - y3
dx31 = x1 - x3
dy31 = y1 - y3
a = dx34*dy31 - dx31*dy34
b = dx12*dy31 - dx31*dy12
d = dx12*dy34 - dx34*dy12
IF (d == 0.) GO TO 1

! D .NE. 0 and the point of intersection of the lines defined by the line
!   segments is P = P1 + (A/D)*(P2-P1) = P3 + (B/D)*(P4-P3).

fn_val = a/d >= 0. .AND. a/d <= 1.  .AND. b/d >= 0. .AND. b/d <= 1.
RETURN

! D .EQ. 0 and thus either the line segments are parallel,
!   or one (or both) of them is a single point.

1 fn_val = a == 0. .AND. b == 0.
RETURN
END FUNCTION intsec



FUNCTION jrand (n) RESULT(fn_val)

! N.B. Arguments IX, IY & IZ have been removed.

INTEGER, INTENT(IN)      :: n
INTEGER                  :: fn_val

!***********************************************************

!                                              From STRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/28/98

!   This function returns a uniformly distributed pseudo-
! random integer in the range 1 to N.


! On input:

!       N = Maximum value to be returned.

! N is not altered by this function.

!       IX,IY,IZ = Integer seeds initialized to values in the range
!                  1 to 30,000 before the first call to JRAND, and not
!                  altered between subsequent calls (unless a sequence of
!                  random numbers is to be repeated by reinitializing the
!                  seeds).

! On output:

!       IX,IY,IZ = Updated integer seeds.

!       JRAND = Random integer in the range 1 to N.

! Reference:  B. A. Wichmann and I. D. Hill, "An Efficient
!             and Portable Pseudo-random Number Generator",
!             Applied Statistics, Vol. 31, No. 2, 1982, pp. 188-190.

! Modules required by JRAND:  None

! Intrinsic functions called by JRAND:  INT, MOD, REAL

!***********************************************************

REAL :: u, x

! Local parameters:

! U = Pseudo-random number uniformly distributed in the interval (0,1).
! X = Pseudo-random number in the range 0 to 3 whose fractional part is U.

ix = MOD(171*ix,30269)
iy = MOD(172*iy,30307)
iz = MOD(170*iz,30323)
x = (REAL(ix)/30269.) + (REAL(iy)/30307.) + (REAL(iz)/30323.)
u = x - INT(x)
fn_val = REAL(n)*u + 1.
RETURN
END FUNCTION jrand



FUNCTION left (x1, y1, x2, y2, x0, y0) RESULT(fn_val)

REAL, INTENT(IN)  :: x1
REAL, INTENT(IN)  :: y1
REAL, INTENT(IN)  :: x2
REAL, INTENT(IN)  :: y2
REAL, INTENT(IN)  :: x0
REAL, INTENT(IN)  :: y0
LOGICAL           :: fn_val

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88

!   This function determines whether node N0 is to the left
! or to the right of the line through N1-N2 as viewed by an
! observer at N1 facing N2.


! On input:

!       X1,Y1 = Coordinates of N1.

!       X2,Y2 = Coordinates of N2.

!       X0,Y0 = Coordinates of N0.

! Input parameters are not altered by this function.

! On output:

!       LEFT = .TRUE. if and only if (X0,Y0) is on or to the
!              left of the directed line N1->N2.

! Modules required by LEFT:  None

!***********************************************************

REAL  :: dx1, dy1, dx2, dy2

! Local parameters:

! DX1,DY1 = X,Y components of the vector N1->N2
! DX2,DY2 = X,Y components of the vector N1->N0

dx1 = x2-x1
dy1 = y2-y1
dx2 = x0-x1
dy2 = y0-y1

! If the sign of the vector cross product of N1->N2 and N1->N0 is positive,
!   then sin(A) > 0, where A is the angle between the vectors, and thus A is
!   in the range (0,180) degrees.

fn_val = dx1*dy2 >= dx2*dy1
RETURN
END FUNCTION left




FUNCTION lstptr (lpl, nb, list, lptr) RESULT(fn_val)

INTEGER, INTENT(IN)  :: lpl
INTEGER, INTENT(IN)  :: nb
INTEGER, INTENT(IN)  :: list(:)
INTEGER, INTENT(IN)  :: lptr(:)
INTEGER              :: fn_val

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88

!   This function returns the index (LIST pointer) of NB in
! the adjacency list for N0, where LPL = LEND(N0).


! On input:

!       LPL = LEND(N0)

!       NB = Index of the node whose pointer is to be returned.
!            NB must be connected to N0.

!       LIST,LPTR = Data structure defining the triangulation
!                   Refer to Subroutine TRMESH.

! Input parameters are not altered by this function.

! On output:

!       LSTPTR = Pointer such that LIST(LSTPTR) = NB or
!                LIST(LSTPTR) = -NB, unless NB is not a
!                neighbor of N0, in which case LSTPTR = LPL.

! Modules required by LSTPTR:  None

!***********************************************************

INTEGER :: lp, nd

lp = lptr(lpl)
1 nd = list(lp)
IF (nd == nb) GO TO 2
lp = lptr(lp)
IF (lp /= lpl) GO TO 1

2 fn_val = lp
RETURN
END FUNCTION lstptr



FUNCTION nbcnt (lpl, lptr) RESULT(fn_val)

INTEGER, INTENT(IN)  :: lpl
INTEGER, INTENT(IN)  :: lptr(:)
INTEGER              :: fn_val

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88

!   This function returns the number of neighbors of a node
! N0 in a triangulation created by Subroutine TRMESH (or TRMSHR).


! On input:

!       LPL = LIST pointer to the last neighbor of N0 --
!             LPL = LEND(N0).

!       LPTR = Array of pointers associated with LIST.

! Input parameters are not altered by this function.

! On output:

!       NBCNT = Number of neighbors of N0.

! Modules required by NBCNT:  None

!***********************************************************

INTEGER :: k, lp

lp = lpl
k = 1

1 lp = lptr(lp)
IF (lp == lpl) GO TO 2
k = k + 1
GO TO 1

2 fn_val = k
RETURN
END FUNCTION nbcnt



FUNCTION nearnd (xp, yp, ist, n, x, y, list, lptr, lend, dsq) RESULT(fn_val)

REAL, INTENT(IN)     :: xp
REAL, INTENT(IN)     :: yp
INTEGER, INTENT(IN)  :: ist
INTEGER, INTENT(IN)  :: n
REAL, INTENT(IN)     :: x(:)
REAL, INTENT(IN)     :: y(:)
INTEGER, INTENT(IN)  :: list(:)
INTEGER, INTENT(IN)  :: lptr(:)
INTEGER, INTENT(IN)  :: lend(:)
REAL, INTENT(OUT)    :: dsq
INTEGER              :: fn_val

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   06/27/98

!   Given a point P in the plane and a Delaunay triangulation created by
! Subroutine TRMESH or TRMSHR, this function returns the index of the nearest
! triangulation node to P.

!   The algorithm consists of implicitly adding P to the triangulation,
! finding the nearest neighbor to P, and implicitly deleting P from the
! triangulation.  Thus, it is based on the fact that, if P is a node in
! a Delaunay triangulation, the nearest node to P is a neighbor of P.


! On input:

!       XP,YP = Cartesian coordinates of the point P to be
!               located relative to the triangulation.

!       IST = Index of a node at which TRFIND begins the search.
!             Search time depends on the proximity of this node to P.

!       N = Number of nodes in the triangulation.  N >= 3.

!       X,Y = Arrays of length N containing the Cartesian
!             coordinates of the nodes.

!       LIST,LPTR,LEND = Data structure defining the triangulation.
!                        Refer to TRMESH.

! Input parameters are not altered by this function.

! On output:

!       NEARND = Nodal index of the nearest node to P, or 0
!                if N < 3 or the triangulation data structure is invalid.

!       DSQ = Squared distance between P and NEARND unless NEARND = 0.

!       Note that the number of candidates for NEARND (neighbors of P)
!       is limited to LMAX defined in the PARAMETER statement below.

! Modules required by NEARND:  JRAND, LEFT, LSTPTR, TRFIND

! Intrinsic function called by NEARND:  ABS

!***********************************************************

INTEGER, PARAMETER :: lmax=25
INTEGER :: i1, i2, i3, l, listp(lmax), lp, lp1, lp2,  &
           lpl, lptrp(lmax), n1, n2, n3, nr, nst
REAL    :: cos1, cos2, ds1, dsr, dx11, dx12, dx21,  &
           dx22, dy11, dy12, dy21, dy22, sin1, sin2

! Store local parameters and test for N invalid.

IF (n < 3) GO TO 7
nst = ist
IF (nst < 1 .OR. nst > n) nst = 1

! Find a triangle (I1,I2,I3) containing P, or the rightmost
!   (I1) and leftmost (I2) visible boundary nodes as viewed from P.

CALL trfind (nst, xp, yp, n, x, y, list, lptr, lend, i1, i2, i3)

! Test for collinear nodes.

IF (i1 == 0) GO TO 7

! Store the linked list of 'neighbors' of P in LISTP and
!   LPTRP.  I1 is the first neighbor, and 0 is stored as
!   the last neighbor if P is not contained in a triangle.
!   L is the length of LISTP and LPTRP, and is limited to LMAX.

IF (i3 /= 0) THEN
  listp(1) = i1
  lptrp(1) = 2
  listp(2) = i2
  lptrp(2) = 3
  listp(3) = i3
  lptrp(3) = 1
  l = 3
ELSE
  n1 = i1
  l = 1
  lp1 = 2
  listp(l) = n1
  lptrp(l) = lp1
  
!   Loop on the ordered sequence of visible boundary nodes N1 from I1 to I2.
  
  1   lpl = lend(n1)
  n1 = -list(lpl)
  l = lp1
  lp1 = l+1
  listp(l) = n1
  lptrp(l) = lp1
  IF (n1 /= i2 .AND. lp1 < lmax) GO TO 1
  l = lp1
  listp(l) = 0
  lptrp(l) = 1
END IF

! Initialize variables for a loop on arcs N1-N2 opposite P in which new
!   'neighbors' are 'swapped' in.  N1 follows N2 as a neighbor of P, and LP1
!   and LP2 are the LISTP indexes of N1 and N2.

lp2 = 1
n2 = i1
lp1 = lptrp(1)
n1 = listp(lp1)

! Begin loop:  find the node N3 opposite N1->N2.

2 lp = lstptr(lend(n1),n2,list,lptr)
IF (list(lp) < 0) GO TO 4
lp = lptr(lp)
n3 = ABS(list(lp))

! Swap test:  Exit the loop if L = LMAX.

IF (l == lmax) GO TO 5
dx11 = x(n1) - x(n3)
dx12 = x(n2) - x(n3)
dx22 = x(n2) - xp
dx21 = x(n1) - xp

dy11 = y(n1) - y(n3)
dy12 = y(n2) - y(n3)
dy22 = y(n2) - yp
dy21 = y(n1) - yp

cos1 = dx11*dx12 + dy11*dy12
cos2 = dx22*dx21 + dy22*dy21
IF (cos1 >= 0. .AND. cos2 >= 0.) GO TO 4
IF (cos1 < 0. .AND. cos2 < 0.) GO TO 3

sin1 = dx11*dy12 - dx12*dy11
sin2 = dx22*dy21 - dx21*dy22
IF (sin1*cos2 + cos1*sin2 >= 0.) GO TO 4

! Swap:  Insert N3 following N2 in the adjacency list for P.
!        The two new arcs opposite P must be tested.

3   l = l+1
lptrp(lp2) = l
listp(l) = n3
lptrp(l) = lp1
lp1 = l
n1 = n3
GO TO 2

! No swap:  Advance to the next arc and test for termination
!           on N1 = I1 (LP1 = 1) or N1 followed by 0.

4   IF (lp1 == 1) GO TO 5
lp2 = lp1
n2 = n1
lp1 = lptrp(lp1)
n1 = listp(lp1)
IF (n1 == 0) GO TO 5
GO TO 2

! Set NR and DSR to the index of the nearest node to P and
!   its squared distance from P, respectively.

5 nr = i1
dsr = (x(nr)-xp)**2 + (y(nr)-yp)**2
DO  lp = 2,l
  n1 = listp(lp)
  IF (n1 == 0) CYCLE
  ds1 = (x(n1)-xp)**2 + (y(n1)-yp)**2
  IF (ds1 < dsr) THEN
    nr = n1
    dsr = ds1
  END IF
END DO
dsq = dsr
fn_val = nr
RETURN

! Invalid input.

7 fn_val = 0
RETURN
END FUNCTION nearnd



SUBROUTINE optim (x, y, na, list, lptr, lend, nit, iwk, ier)

REAL, INTENT(IN)         :: x(:)
REAL, INTENT(IN)         :: y(:)
INTEGER, INTENT(IN)      :: na
INTEGER, INTENT(IN OUT)  :: list(:)
INTEGER, INTENT(IN OUT)  :: lptr(:)
INTEGER, INTENT(IN OUT)  :: lend(:)
INTEGER, INTENT(IN OUT)  :: nit
INTEGER, INTENT(IN OUT)  :: iwk(:,:)
INTEGER, INTENT(OUT)     :: ier

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   06/27/98

!   Given a set of NA triangulation arcs, this subroutine
! optimizes the portion of the triangulation consisting of
! the quadrilaterals (pairs of adjacent triangles) which
! have the arcs as diagonals by applying the circumcircle
! test and appropriate swaps to the arcs.

!   An iteration consists of applying the swap test and
! swaps to all NA arcs in the order in which they are
! stored.  The iteration is repeated until no swap occurs
! or NIT iterations have been performed.  The bound on the
! number of iterations may be necessary to prevent an
! infinite loop caused by cycling (reversing the effect of a
! previous swap) due to floating point inaccuracy when four
! or more nodes are nearly cocircular.


! On input:

!       X,Y = Arrays containing the nodal coordinates.

!       NA = Number of arcs in the set.  NA >= 0.

! The above parameters are not altered by this routine.

!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to Subroutine TRMESH.

!       NIT = Maximum number of iterations to be performed.
!             A reasonable value is 3*NA.  NIT >= 1.

!       IWK = Integer array dimensioned 2 by NA containing the nodal indexes
!             of the arc endpoints (pairs of endpoints are stored in columns).

! On output:

!       LIST,LPTR,LEND = Updated triangulation data structure reflecting the
!                        swaps.

!       NIT = Number of iterations performed.

!       IWK = Endpoint indexes of the new set of arcs reflecting the swaps.

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if a swap occurred on the last of
!                     MAXIT iterations, where MAXIT is the
!                     value of NIT on input.  The new set
!                     of arcs in not necessarily optimal in this case.
!             IER = 2 if NA < 0 or NIT < 1 on input.
!             IER = 3 if IWK(2,I) is not a neighbor of
!                     IWK(1,I) for some I in the range 1
!                     to NA.  A swap may have occurred in this case.
!             IER = 4 if a zero pointer was returned by Subroutine SWAP.

! Modules required by OPTIM:  LSTPTR, SWAP, SWPTST

! Intrinsic function called by OPTIM:  ABS

!***********************************************************

INTEGER :: i, io1, io2, iter, lp, lp21, lpl, lpp, maxit, n1, n2, nna
LOGICAL :: swp

! Local parameters:

! I =       Column index for IWK
! IO1,IO2 = Nodal indexes of the endpoints of an arc in IWK
! ITER =    Iteration count
! LP =      LIST pointer
! LP21 =    Parameter returned by SWAP (not used)
! LPL =     Pointer to the last neighbor of IO1
! LPP =     Pointer to the node preceding IO2 as a neighbor of IO1
! MAXIT =   Input value of NIT
! N1,N2 =   Nodes opposite IO1->IO2 and IO2->IO1, respectively
! NNA =     Local copy of NA
! SWP =     Flag set to TRUE iff a swap occurs in the optimization loop

nna = na
maxit = nit
IF (nna < 0 .OR. maxit < 1) GO TO 7

! Initialize iteration count ITER and test for NA = 0.

iter = 0
IF (nna == 0) GO TO 5

! Top of loop --
!   SWP = TRUE iff a swap occurred in the current iteration.

1 IF (iter == maxit) GO TO 6
iter = iter + 1
swp = .false.

!   Inner loop on arcs IO1-IO2 --

DO  i = 1,nna
  io1 = iwk(1,i)
  io2 = iwk(2,i)
  
!   Set N1 and N2 to the nodes opposite IO1->IO2 and
!     IO2->IO1, respectively.  Determine the following:
  
!     LPL = pointer to the last neighbor of IO1,
!     LP = pointer to IO2 as a neighbor of IO1, and
!     LPP = pointer to the node N2 preceding IO2.
  
  lpl = lend(io1)
  lpp = lpl
  lp = lptr(lpp)
  2   IF (list(lp) == io2) GO TO 3
  lpp = lp
  lp = lptr(lpp)
  IF (lp /= lpl) GO TO 2
  
!   IO2 should be the last neighbor of IO1.  Test for no
!     arc and bypass the swap test if IO1 is a boundary node.
  
  IF (ABS(list(lp)) /= io2) GO TO 8
  IF (list(lp) < 0) CYCLE
  
!   Store N1 and N2, or bypass the swap test if IO1 is a
!     boundary node and IO2 is its first neighbor.
  
  3   n2 = list(lpp)
  IF (n2 < 0) CYCLE
  lp = lptr(lp)
  n1 = ABS(list(lp))
  
!   Test IO1-IO2 for a swap, and update IWK if necessary.
  
  IF ( .NOT. swptst(n1,n2,io1,io2,x,y) ) CYCLE
  CALL swap (n1,n2,io1,io2, list,lptr,lend, lp21)
  IF (lp21 == 0) GO TO 9
  swp = .true.
  iwk(1,i) = n1
  iwk(2,i) = n2
END DO
IF (swp) GO TO 1

! Successful termination.

5 nit = iter
ier = 0
RETURN

! MAXIT iterations performed without convergence.

6 nit = maxit
ier = 1
RETURN

! Invalid input parameter.

7 nit = 0
ier = 2
RETURN

! IO2 is not a neighbor of IO1.

8 nit = iter
ier = 3
RETURN

! Zero pointer returned by SWAP.

9 nit = iter
ier = 4
RETURN
END SUBROUTINE optim



FUNCTION store (x) RESULT(fn_val)

REAL, INTENT(IN)  :: x
REAL              :: fn_val

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   03/18/90

!   This function forces its argument X to be stored in a memory location,
! thus providing a means of determining floating point number characteristics
! (such as the machine precision) when it is necessary to avoid computation
! in high precision registers.

! On input:

!       X = Value to be stored.

! X is not altered by this function.

! On output:

!       STORE = Value of X after it has been stored and
!               possibly truncated or rounded to the single
!               precision word length.

! Modules required by STORE:  None

!***********************************************************

stcom = x
fn_val = stcom
RETURN
END FUNCTION store



SUBROUTINE swap (in1, in2, io1, io2, list, lptr, lend, lp21)

INTEGER, INTENT(IN)      :: in1
INTEGER, INTENT(IN)      :: in2
INTEGER, INTENT(IN)      :: io1
INTEGER, INTENT(IN)      :: io2
INTEGER, INTENT(IN OUT)  :: list(:)
INTEGER, INTENT(IN OUT)  :: lptr(:)
INTEGER, INTENT(IN OUT)  :: lend(:)
INTEGER, INTENT(OUT)     :: lp21

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   06/22/98

!   Given a triangulation of a set of points on the unit sphere,
! this subroutine replaces a diagonal arc in a strictly convex
! quadrilateral (defined by a pair of adjacent triangles) with the
! other diagonal.  Equivalently, a pair of adjacent triangles is
! replaced by another pair having the same union.

! On input:

!       IN1,IN2,IO1,IO2 = Nodal indexes of the vertices of
!                         the quadrilateral.  IO1-IO2 is re-
!                         placed by IN1-IN2.  (IO1,IO2,IN1)
!                         and (IO2,IO1,IN2) must be trian-
!                         gles on input.

! The above parameters are not altered by this routine.

!       LIST,LPTR,LEND = Data structure defining the triangulation.
!                        Refer to Subroutine TRMESH.

! On output:

!       LIST,LPTR,LEND = Data structure updated with the swap --
!                        triangles (IO1,IO2,IN1) and (IO2,IO1,IN2) are
!                        replaced by (IN1,IN2,IO2) and (IN2,IN1,IO1)
!                        unless LP21 = 0.

!       LP21 = Index of IN1 as a neighbor of IN2 after the
!              swap is performed unless IN1 and IN2 are
!              adjacent on input, in which case LP21 = 0.

! Module required by SWAP:  LSTPTR

! Intrinsic function called by SWAP:  ABS

!***********************************************************

INTEGER :: lp, lph, lpsav

! Local parameters:

! LP,LPH,LPSAV = LIST pointers


! Test for IN1 and IN2 adjacent.

lp = lstptr(lend(in1), in2, list, lptr)
IF (ABS(list(lp)) == in2) THEN
  lp21 = 0
  RETURN
END IF

! Delete IO2 as a neighbor of IO1.

lp = lstptr(lend(io1), in2, list, lptr)
lph = lptr(lp)
lptr(lp) = lptr(lph)

! If IO2 is the last neighbor of IO1, make IN2 the last neighbor.

IF (lend(io1) == lph) lend(io1) = lp

! Insert IN2 as a neighbor of IN1 following IO1 using the hole created above.

lp = lstptr(lend(in1), io1, list, lptr)
lpsav = lptr(lp)
lptr(lp) = lph
list(lph) = in2
lptr(lph) = lpsav

! Delete IO1 as a neighbor of IO2.

lp = lstptr(lend(io2), in1, list, lptr)
lph = lptr(lp)
lptr(lp) = lptr(lph)

! If IO1 is the last neighbor of IO2, make IN1 the last neighbor.

IF (lend(io2) == lph) lend(io2) = lp

! Insert IN1 as a neighbor of IN2 following IO2.

lp = lstptr(lend(in2),io2,list,lptr)
lpsav = lptr(lp)
lptr(lp) = lph
list(lph) = in1
lptr(lph) = lpsav
lp21 = lph
RETURN
END SUBROUTINE swap



FUNCTION swptst (in1, in2, io1, io2, x, y) RESULT(fn_val)

INTEGER, INTENT(IN)  :: in1
INTEGER, INTENT(IN)  :: in2
INTEGER, INTENT(IN)  :: io1
INTEGER, INTENT(IN)  :: io2
REAL, INTENT(IN)     :: x(:)
REAL, INTENT(IN)     :: y(:)
LOGICAL              :: fn_val

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88

!   This function applies the circumcircle test to a quadri-
! lateral defined by a pair of adjacent triangles.  The
! diagonal arc (shared triangle side) should be swapped for
! the other diagonl if and only if the fourth vertex is
! strictly interior to the circumcircle of one of the
! triangles (the decision is independent of the choice of
! triangle).  Equivalently, the diagonal is chosen to maxi-
! mize the smallest of the six interior angles over the two
! pairs of possible triangles (the decision is for no swap
! if the quadrilateral is not strictly convex).

!   When the four vertices are nearly cocircular (the
! neutral case), the preferred decision is no swap -- in
! order to avoid unnecessary swaps and, more important, to
! avoid cycling in Subroutine OPTIM which is called by
! DELNOD and EDGE.  Thus, a tolerance SWTOL (stored in
! SWPCOM by TRMESH or TRMSHR) is used to define 'nearness'
! to the neutral case.


! On input:

!       IN1,IN2,IO1,IO2 = Nodal indexes of the vertices of
!                         the quadrilateral.  IO1-IO2 is the
!                         triangulation arc (shared triangle
!                         side) to be replaced by IN1-IN2 if
!                         the decision is to swap.  The
!                         triples (IO1,IO2,IN1) and (IO2,
!                         IO1,IN2) must define triangles (be
!                         in counterclockwise order) on input.

!       X,Y = Arrays containing the nodal coordinates.

! Input parameters are not altered by this routine.

! On output:

!       SWPTST = .TRUE. if and only if the arc connecting
!                IO1 and IO2 is to be replaced.

! Modules required by SWPTST:  None

!***********************************************************

REAL :: dx11, dx12, dx22, dx21, dy11, dy12, dy22, dy21,  &
        sin1, sin2, cos1, cos2, sin12

! Local parameters:

! DX11,DY11 = X,Y components of the vector IN1->IO1
! DX12,DY12 = X,Y components of the vector IN1->IO2
! DX22,DY22 = X,Y components of the vector IN2->IO2
! DX21,DY21 = X,Y components of the vector IN2->IO1
! SIN1 =      Cross product of the vectors IN1->IO1 and
!               IN1->IO2 -- proportional to sin(T1), where
!               T1 is the angle at IN1 formed by the vectors
! COS1 =      Inner product of the vectors IN1->IO1 and
!               IN1->IO2 -- proportional to cos(T1)
! SIN2 =      Cross product of the vectors IN2->IO2 and
!               IN2->IO1 -- proportional to sin(T2), where
!               T2 is the angle at IN2 formed by the vectors
! COS2 =      Inner product of the vectors IN2->IO2 and
!               IN2->IO1 -- proportional to cos(T2)
! SIN12 =     SIN1*COS2 + COS1*SIN2 -- proportional to
!               sin(T1+T2)


! Compute the vectors containing the angles T1 and T2.

dx11 = x(io1) - x(in1)
dx12 = x(io2) - x(in1)
dx22 = x(io2) - x(in2)
dx21 = x(io1) - x(in2)

dy11 = y(io1) - y(in1)
dy12 = y(io2) - y(in1)
dy22 = y(io2) - y(in2)
dy21 = y(io1) - y(in2)

! Compute inner products.

cos1 = dx11*dx12 + dy11*dy12
cos2 = dx22*dx21 + dy22*dy21

! The diagonals should be swapped iff (T1+T2) > 180 degrees.
!   The following two tests ensure numerical stability:
!   the decision must be FALSE when both angles are close to 0,
!   and TRUE when both angles are close to 180 degrees.

IF (cos1 >= 0. .AND. cos2 >= 0.) GO TO 2
IF (cos1 < 0. .AND. cos2 < 0.) GO TO 1

! Compute vector cross products (Z-components).

sin1 = dx11*dy12 - dx12*dy11
sin2 = dx22*dy21 - dx21*dy22
sin12 = sin1*cos2 + cos1*sin2
IF (sin12 >= -swtol) GO TO 2

! Swap.

1 fn_val = .true.
RETURN

! No swap.

2 fn_val = .false.
RETURN
END FUNCTION swptst



SUBROUTINE trfind (nst, px, py, n, x, y, list, lptr, lend, i1, i2, i3)

INTEGER, INTENT(IN)   :: nst
REAL, INTENT(IN)      :: px
REAL, INTENT(IN)      :: py
INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN)      :: x(:)
REAL, INTENT(IN)      :: y(:)
INTEGER, INTENT(IN)   :: list(:)
INTEGER, INTENT(IN)   :: lptr(:)
INTEGER, INTENT(IN)   :: lend(:)
INTEGER, INTENT(OUT)  :: i1
INTEGER, INTENT(OUT)  :: i2
INTEGER, INTENT(OUT)  :: i3

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/28/98

!   This subroutine locates a point P relative to a triangu-
! lation created by Subroutine TRMESH or TRMSHR.  If P is
! contained in a triangle, the three vertex indexes are
! returned.  Otherwise, the indexes of the rightmost and
! leftmost visible boundary nodes are returned.


! On input:

!       NST = Index of a node at which TRFIND begins the
!             search.  Search time depends on the proximity
!             of this node to P.

!       PX,PY = X and y coordinates of the point P to be located.

!       N = Number of nodes in the triangulation.  N >= 3.

!       X,Y = Arrays of length N containing the coordinates
!             of the nodes in the triangulation.

!       LIST,LPTR,LEND = Data structure defining the triangulation.
!                        Refer to Subroutine TRMESH.

! Input parameters are not altered by this routine.

! On output:

!       I1,I2,I3 = Nodal indexes, in counterclockwise order,
!                  of the vertices of a triangle containing
!                  P if P is contained in a triangle.  If P
!                  is not in the convex hull of the nodes,
!                  I1 indexes the rightmost visible boundary
!                  node, I2 indexes the leftmost visible
!                  boundary node, and I3 = 0.  Rightmost and
!                  leftmost are defined from the perspective
!                  of P, and a pair of points are visible
!                  from each other if and only if the line
!                  segment joining them intersects no trian-
!                  gulation arc.  If P and all of the nodes
!                  lie on a common line, then I1 = I2 = I3 =
!                  0 on output.

! Modules required by TRFIND:  JRAND, LEFT, LSTPTR, STORE

! Intrinsic function called by TRFIND:  ABS

!***********************************************************

INTEGER :: lp, n0, n1, n1s, n2, n2s, n3, n4, nb, nf, nl, np, npp
REAL    :: b1, b2, xp, yp

! Local parameters:

! B1,B2 =    Unnormalized barycentric coordinates of P with
!              respect to (N1,N2,N3)
! IX,IY,IZ = Integer seeds for JRAND
! LP =       LIST pointer
! N0,N1,N2 = Nodes in counterclockwise order defining a
!              cone (with vertex N0) containing P
! N1S,N2S =  Saved values of N1 and N2
! N3,N4 =    Nodes opposite N1->N2 and N2->N1, respectively
! NB =       Index of a boundary node -- first neighbor of
!              NF or last neighbor of NL in the boundary traversal loops
! NF,NL =    First and last neighbors of N0, or first
!              (rightmost) and last (leftmost) nodes
!              visible from P when P is exterior to the triangulation
! NP,NPP =   Indexes of boundary nodes used in the boundary traversal loops
! XA,XB,XC = Dummy arguments for FRWRD
! YA,YB,YC = Dummy arguments for FRWRD
! XP,YP =    Local variables containing the components of P

! Statement function:

! FRWRD = TRUE iff C is forward of A->B
!              iff <A->B,A->C> >= 0.

! frwrd(xa,ya,xb,yb,xc,yc) = (xb-xa)*(xc-xa) + (yb-ya)*(yc-ya) >= 0.

! Initialize variables.

xp = px
yp = py
n0 = nst
IF (n0 < 1 .OR. n0 > n) n0 = jrand(n)

! Set NF and NL to the first and last neighbors of N0, and initialize N1 = NF.

1 lp = lend(n0)
nl = list(lp)
lp = lptr(lp)
nf = list(lp)
n1 = nf

! Find a pair of adjacent neighbors N1,N2 of N0 that define
!   a wedge containing P:  P LEFT N0->N1 and P RIGHT N0->N2.

IF (nl > 0) GO TO 2

!   N0 is a boundary node.  Test for P exterior.

nl = -nl
IF ( .NOT. left(x(n0),y(n0),x(nf),y(nf),xp,yp) ) THEN
  nl = n0
  GO TO 9
END IF
IF ( .NOT. left(x(nl),y(nl),x(n0),y(n0),xp,yp) ) THEN
  nb = nf
  nf = n0
  np = nl
  npp = n0
  GO TO 11
END IF
GO TO 3

!   N0 is an interior node.  Find N1.

2 IF ( left(x(n0), y(n0), x(n1), y(n1), xp, yp) ) GO TO 3
lp = lptr(lp)
n1 = list(lp)
IF (n1 == nl) GO TO 6
GO TO 2

!   P is to the left of edge N0->N1.  Initialize N2 to the
!     next neighbor of N0.

3 lp = lptr(lp)
n2 = ABS(list(lp))
IF ( .NOT. left(x(n0),y(n0),x(n2),y(n2),xp,yp) ) GO TO 7
n1 = n2
IF (n1 /= nl) GO TO 3
IF ( .NOT. left(x(n0),y(n0),x(nf),y(nf),xp,yp) ) GO TO 6
IF (xp == x(n0) .AND. yp == y(n0)) GO TO 5

!   P is left of or on edges N0->NB for all neighbors NB of N0.
!   All points are collinear iff P is left of NB->N0 for
!     all neighbors NB of N0.  Search the neighbors of N0.
!     NOTE -- N1 = NL and LP points to NL.

4 IF ( .NOT. left(x(n1), y(n1), x(n0), y(n0), xp, yp) ) GO TO 5
lp = lptr(lp)
n1 = ABS(list(lp))
IF (n1 == nl) GO TO 17
GO TO 4

!   P is to the right of N1->N0, or P=N0.  Set N0 to N1 and start over.

5 n0 = n1
GO TO 1

!   P is between edges N0->N1 and N0->NF.

6 n2 = nf

! P is contained in the wedge defined by line segments
!   N0->N1 and N0->N2, where N1 is adjacent to N2.  Set
!   N3 to the node opposite N1->N2, and save N1 and N2 to
!   test for cycling.

7 n3 = n0
n1s = n1
n2s = n2

! Top of edge hopping loop.  Test for termination.

8 IF ( left(x(n1), y(n1), x(n2), y(n2), xp, yp) ) THEN
  
!   P LEFT N1->N2 and hence P is in (N1,N2,N3) unless an
!     error resulted from floating point inaccuracy and
!     collinearity.  Compute the unnormalized barycentric
!     coordinates of P with respect to (N1,N2,N3).
  
  b1 = (x(n3)-x(n2))*(yp-y(n2)) - (xp-x(n2))*(y(n3)-y(n2))
  b2 = (x(n1)-x(n3))*(yp-y(n3)) - (xp-x(n3))*(y(n1)-y(n3))
  IF (store(b1+1.) >= 1.  .AND. store(b2+1.) >= 1.) GO TO 16
  
!   Restart with N0 randomly selected.
  
  n0 = jrand(n)
  GO TO 1
END IF

!   Set N4 to the neighbor of N2 which follows N1 (node
!     opposite N2->N1) unless N1->N2 is a boundary edge.

lp = lstptr(lend(n2), n1, list, lptr)
IF (list(lp) < 0) THEN
  nf = n2
  nl = n1
  GO TO 9
END IF
lp = lptr(lp)
n4 = ABS(list(lp))

!   Select the new edge N1->N2 which intersects the line
!     segment N0-P, and set N3 to the node opposite N1->N2.

IF ( left(x(n0), y(n0), x(n4), y(n4), xp, yp) ) THEN
  n3 = n1
  n1 = n4
  n2s = n2
  IF (n1 /= n1s .AND. n1 /= n0) GO TO 8
ELSE
  n3 = n2
  n2 = n4
  n1s = n1
  IF (n2 /= n2s .AND. n2 /= n0) GO TO 8
END IF

!   The starting node N0 or edge N1-N2 was encountered
!     again, implying a cycle (infinite loop).  Restart
!     with N0 randomly selected.

n0 = jrand(n)
GO TO 1

! Boundary traversal loops.  NL->NF is a boundary edge and
!   P RIGHT NL->NF.  Save NL and NF.

9 np = nl
npp = nf

! Find the first (rightmost) visible boundary node NF.  NB
!   is set to the first neighbor of NF, and NP is the last neighbor.

10 lp = lend(nf)
lp = lptr(lp)
nb = list(lp)
IF ( .NOT. left(x(nf),y(nf),x(nb),y(nb),xp,yp) ) GO TO 12

!   P LEFT NF->NB and thus NB is not visible unless an error
!     resulted from floating point inaccuracy and collinear-
!     ity of the 4 points NP, NF, NB, and P.

11 IF ( frwrd(x(nf),y(nf),x(np),y(np),xp,yp) .OR. &
      frwrd(x(nf),y(nf),x(np),y(np),x(nb),y(nb)) ) THEN
  i1 = nf
  GO TO 13
END IF

!   Bottom of loop.

12 np = nf
nf = nb
GO TO 10

! Find the last (leftmost) visible boundary node NL.  NB
!   is set to the last neighbor of NL, and NPP is the first neighbor.

13 lp = lend(nl)
nb = -list(lp)
IF ( .NOT. left(x(nb),y(nb),x(nl),y(nl),xp,yp) ) GO TO 14

!   P LEFT NB->NL and thus NB is not visible unless an error
!     resulted from floating point inaccuracy and collinear-
!     ity of the 4 points P, NB, NL, and NPP.

IF ( frwrd(x(nl),y(nl),x(npp),y(npp),xp,yp) .OR. &
    frwrd(x(nl),y(nl),x(npp),y(npp),x(nb),y(nb)) ) GO TO 15

!   Bottom of loop.

14 npp = nl
nl = nb
GO TO 13

! NL is the leftmost visible boundary node.

15 i2 = nl
i3 = 0
RETURN

! P is in the triangle (N1,N2,N3).

16 i1 = n1
i2 = n2
i3 = n3
RETURN

! All points are collinear.

17 i1 = 0
i2 = 0
i3 = 0
RETURN
END SUBROUTINE trfind



SUBROUTINE trlist (ncc, lcc, n, list, lptr, lend, nrow, nt, ltri, lct, ier)

INTEGER, INTENT(IN)   :: ncc
INTEGER, INTENT(IN)   :: lcc(:)
INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: list(:)
INTEGER, INTENT(IN)   :: lptr(:)
INTEGER, INTENT(IN)   :: lend(:)
INTEGER, INTENT(IN)   :: nrow
INTEGER, INTENT(OUT)  :: nt
INTEGER, INTENT(OUT)  :: ltri(:,:)   ! ltri(nrow,*)
INTEGER, INTENT(OUT)  :: lct(:)
INTEGER, INTENT(OUT)  :: ier

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   03/22/97

!   This subroutine converts a triangulation data structure
! from the linked list created by Subroutine TRMESH or
! TRMSHR to a triangle list.

! On input:

!       NCC = Number of constraints.  NCC >= 0.

!       LCC = List of constraint curve starting indexes (or
!             dummy array of length 1 if NCC = 0).  Refer to
!             Subroutine ADDCST.

!       N = Number of nodes in the triangulation.  N >= 3.

!       LIST,LPTR,LEND = Linked list data structure defin-
!                        ing the triangulation.  Refer to
!                        Subroutine TRMESH.

!       NROW = Number of rows (entries per triangle) reserved for the
!              triangle list LTRI.  The value must be 6 if only the vertex
!              indexes and neighboring triangle indexes are to be stored,
!              or 9 if arc indexes are also to be assigned and stored.
!              Refer to LTRI.

! The above parameters are not altered by this routine.

!       LTRI = Integer array of length at least NROW*NT,
!              where NT is at most 2N-5.  (A sufficient
!              length is 12N if NROW=6 or 18N if NROW=9.)

!       LCT = Integer array of length NCC or dummy array of
!             length 1 if NCC = 0.

! On output:

!       NT = Number of triangles in the triangulation unless
!            IER .NE. 0, in which case NT = 0.  NT = 2N - NB
!            - 2, where NB is the number of boundary nodes.

!       LTRI = NROW by NT array whose J-th column contains
!              the vertex nodal indexes (first three rows),
!              neighboring triangle indexes (second three
!              rows), and, if NROW = 9, arc indexes (last
!              three rows) associated with triangle J for
!              J = 1,...,NT.  The vertices are ordered
!              counterclockwise with the first vertex taken
!              to be the one with smallest index.  Thus,
!              LTRI(2,J) and LTRI(3,J) are larger than
!              LTRI(1,J) and index adjacent neighbors of
!              node LTRI(1,J).  For I = 1,2,3, LTRI(I+3,J)
!              and LTRI(I+6,J) index the triangle and arc,
!              respectively, which are opposite (not shared
!              by) node LTRI(I,J), with LTRI(I+3,J) = 0 if
!              LTRI(I+6,J) indexes a boundary arc.  Vertex
!              indexes range from 1 to N, triangle indexes
!              from 0 to NT, and, if included, arc indexes
!              from 1 to NA = NT+N-1.  The triangles are or-
!              dered on first (smallest) vertex indexes,
!              except that the sets of constraint triangles
!              (triangles contained in the closure of a con-
!              straint region) follow the non-constraint triangles.

!       LCT = Array of length NCC containing the triangle
!             index of the first triangle of constraint J in
!             LCT(J).  Thus, the number of non-constraint
!             triangles is LCT(1)-1, and constraint J con-
!             tains LCT(J+1)-LCT(J) triangles, where
!             LCT(NCC+1) = NT+1.

!       IER = Error indicator.
!             IER = 0 if no errors were encountered.
!             IER = 1 if NCC, N, NROW, or an LCC entry is
!                     outside its valid range on input.
!             IER = 2 if the triangulation data structure
!                     (LIST,LPTR,LEND) is invalid.  Note,
!                     however, that these arrays are not
!                     completely tested for validity.

! Modules required by TRLIST:  None

! Intrinsic function called by TRLIST:  ABS

!***********************************************************

INTEGER :: i, i1, i2, i3, isv, j, jlast, ka, kn, kt, l,  &
           lcc1, lp, lp2, lpl, lpln1, n1, n1st, n2, n3, nm2, nn
LOGICAL :: arcs, cstri, pass2

! Test for invalid input parameters and store the index
!   LCC1 of the first constraint node (if any).

nn = n
IF (ncc < 0 .OR. (nrow /= 6  .AND. nrow /= 9)) GO TO 12
lcc1 = nn+1
IF (ncc == 0) THEN
  IF (nn < 3) GO TO 12
ELSE
  DO  i = ncc,1,-1
    IF (lcc1-lcc(i) < 3) GO TO 12
    lcc1 = lcc(i)
  END DO
  IF (lcc1 < 1) GO TO 12
END IF

! Initialize parameters for loop on triangles KT = (N1,N2,
!   N3), where N1 < N2 and N1 < N3.  This requires two
!   passes through the nodes with all non-constraint
!   triangles stored on the first pass, and the constraint
!   triangles stored on the second.

!   ARCS = TRUE iff arc indexes are to be stored.
!   KA,KT = Numbers of currently stored arcs and triangles.
!   N1ST = Starting index for the loop on nodes (N1ST = 1 on
!            pass 1, and N1ST = LCC1 on pass 2).
!   NM2 = Upper bound on candidates for N1.
!   PASS2 = TRUE iff constraint triangles are to be stored.

arcs = nrow == 9
ka = 0
kt = 0
n1st = 1
nm2 = nn-2
pass2 = .false.

! Loop on nodes N1:  J = constraint containing N1,
!                    JLAST = last node in constraint J.

2 j = 0
jlast = lcc1 - 1
DO  n1 = n1st,nm2
  IF (n1 > jlast) THEN
    
! N1 is the first node in constraint J+1.  Update J and
!   JLAST, and store the first constraint triangle index
!   if in pass 2.
    
    j = j + 1
    IF (j < ncc) THEN
      jlast = lcc(j+1) - 1
    ELSE
      jlast = nn
    END IF
    IF (pass2) lct(j) = kt + 1
  END IF
  
! Loop on pairs of adjacent neighbors (N2,N3).  LPLN1 points
!   to the last neighbor of N1, and LP2 points to N2.
  
  lpln1 = lend(n1)
  lp2 = lpln1
  3     lp2 = lptr(lp2)
  n2 = list(lp2)
  lp = lptr(lp2)
  n3 = ABS(list(lp))
  IF (n2 < n1 .OR. n3 < n1) GO TO 10
  
! (N1,N2,N3) is a constraint triangle iff the three nodes
!   are in the same constraint and N2 < N3.  Bypass con-
!   straint triangles on pass 1 and non-constraint triangles
!   on pass 2.
  
  cstri = n1 >= lcc1 .AND. n2 < n3  .AND. n3 <= jlast
  IF ((cstri .AND. .NOT. pass2)  .OR. (.NOT. cstri .AND. pass2)) GO TO 10
  
! Add a new triangle KT = (N1,N2,N3).
  
  kt = kt + 1
  ltri(1,kt) = n1
  ltri(2,kt) = n2
  ltri(3,kt) = n3
  
! Loop on triangle sides (I1,I2) with neighboring triangles
!   KN = (I1,I2,I3).
  
  DO  i = 1,3
    IF (i == 1) THEN
      i1 = n3
      i2 = n2
    ELSE IF (i == 2) THEN
      i1 = n1
      i2 = n3
    ELSE
      i1 = n2
      i2 = n1
    END IF
    
! Set I3 to the neighbor of I1 which follows I2 unless
!   I2->I1 is a boundary arc.
    
    lpl = lend(i1)
    lp = lptr(lpl)
    4 IF (list(lp) == i2) GO TO 5
    lp = lptr(lp)
    IF (lp /= lpl) GO TO 4
    
!   I2 is the last neighbor of I1 unless the data structure
!     is invalid.  Bypass the search for a neighboring
!     triangle if I2->I1 is a boundary arc.
    
    IF (ABS(list(lp)) /= i2) GO TO 13
    kn = 0
    IF (list(lp) < 0) GO TO 8
    
!   I2->I1 is not a boundary arc, and LP points to I2 as a neighbor of I1.
    
    5 lp = lptr(lp)
    i3 = ABS(list(lp))
    
! Find L such that LTRI(L,KN) = I3 (not used if KN > KT),
!   and permute the vertex indexes of KN so that I1 is
!   smallest.
    
    IF (i1 < i2 .AND. i1 < i3) THEN
      l = 3
    ELSE IF (i2 < i3) THEN
      l = 2
      isv = i1
      i1 = i2
      i2 = i3
      i3 = isv
    ELSE
      l = 1
      isv = i1
      i1 = i3
      i3 = i2
      i2 = isv
    END IF
    
! Test for KN > KT (triangle index not yet assigned).
    
    IF (i1 > n1 .AND. .NOT. pass2) CYCLE
    
! Find KN, if it exists, by searching the triangle list in
!   reverse order.
    
    DO  kn = kt-1,1,-1
      IF (ltri(1,kn) == i1 .AND. ltri(2,kn) ==  &
          i2 .AND. ltri(3,kn) == i3) GO TO 7
    END DO
    CYCLE
    
! Store KT as a neighbor of KN.
    
    7 ltri(l+3,kn) = kt
    
! Store KN as a neighbor of KT, and add a new arc KA.
    
    8 ltri(i+3,kt) = kn
    IF (arcs) THEN
      ka = ka + 1
      ltri(i+6,kt) = ka
      IF (kn /= 0) ltri(l+6,kn) = ka
    END IF
  END DO
  
! Bottom of loop on triangles.
  
  10 IF (lp2 /= lpln1) GO TO 3
END DO

! Bottom of loop on nodes.

IF (.NOT. pass2 .AND. ncc > 0) THEN
  pass2 = .true.
  n1st = lcc1
  GO TO 2
END IF

! No errors encountered.

nt = kt
ier = 0
RETURN

! Invalid input parameter.

12 nt = 0
ier = 1
RETURN

! Invalid triangulation data structure:  I1 is a neighbor of
!   I2, but I2 is not a neighbor of I1.

13 nt = 0
ier = 2
RETURN
END SUBROUTINE trlist



SUBROUTINE trlprt (ncc, lct, n, x, y, nrow, nt, ltri, lout, prntx)

INTEGER, INTENT(IN)      :: ncc
INTEGER, INTENT(IN OUT)  :: lct(:)
INTEGER, INTENT(IN)      :: n
REAL, INTENT(IN OUT)     :: x(:)
REAL, INTENT(IN OUT)     :: y(:)
INTEGER, INTENT(IN)      :: nrow
INTEGER, INTENT(IN)      :: nt
INTEGER, INTENT(IN OUT)  :: ltri(:,:)   ! ltri(nrow,nt)
INTEGER, INTENT(IN)      :: lout
LOGICAL, INTENT(IN)      :: prntx

!***********************************************************

!                                               From TRLPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/02/98

!   Given a triangulation of a set of points in the plane,
! this subroutine prints the triangle list created by
! Subroutine TRLIST and, optionally, the nodal coordinates
! on logical unit LOUT.  The numbers of boundary nodes,
! triangles, and arcs, and the constraint region triangle
! indexes, if any, are also printed.

!   All parameters other than LOUT and PRNTX should be
! unaltered from their values on output from TRLIST.


! On input:

!       NCC = Number of constraints.

!       LCT = List of constraint triangle starting indexes
!             (or dummy array of length 1 if NCC = 0).

!       N = Number of nodes in the triangulation.
!           3 <= N <= 9999.

!       X,Y = Arrays of length N containing the coordinates
!             of the nodes in the triangulation -- not used
!             unless PRNTX = TRUE.

!       NROW = Number of rows (entries per triangle) re-
!              served for the triangle list LTRI.  The value
!              must be 6 if only the vertex indexes and
!              neighboring triangle indexes are stored, or 9
!              if arc indexes are also stored.

!       NT = Number of triangles in the triangulation.
!            1 <= NT <= 9999.

!       LTRI = NROW by NT array whose J-th column contains
!              the vertex nodal indexes (first three rows),
!              neighboring triangle indexes (second three
!              rows), and, if NROW = 9, arc indexes (last
!              three rows) associated with triangle J for
!              J = 1,...,NT.

!       LOUT = Logical unit number for output.  0 <= LOUT
!              <= 99.  Output is printed on unit 6 if LOUT
!              is outside its valid range on input.

!       PRNTX = Logical variable with value TRUE if and only
!               if X and Y are to be printed (to 6 decimal places).

! None of the parameters are altered by this routine.

! Modules required by TRLPRT:  None

!***********************************************************

INTEGER :: i, k, lun, na, nb, nl
INTEGER, PARAMETER  :: nmax = 9999, nlmax = 60

! Local parameters:

!   I = DO-loop, nodal index, and row index for LTRI
!   K = DO-loop and triangle index
!   LUN = Logical unit number for output
!   NA = Number of triangulation arcs
!   NB = Number of boundary nodes
!   NL = Number of lines printed on the current page
!   NLMAX = Maximum number of print lines per page
!   NMAX = Maximum value of N and NT (4-digit format)

lun = lout
IF (lun < 0 .OR. lun > 99) lun = 6

! Print a heading and test for invalid input.

WRITE (lun,100)
nl = 1
IF (n < 3 .OR. n > nmax  .OR. (nrow /= 6 .AND. nrow /= 9) .OR. &
      nt < 1 .OR. nt > nmax) THEN
  
! Print an error message and bypass the loops.
  
  WRITE (lun,110) n, nrow, nt
  GO TO 3
END IF
IF (prntx) THEN
  
! Print X and Y.
  
  WRITE (lun,101)
  nl = 6
  DO  i = 1,n
    IF (nl >= nlmax) THEN
      WRITE (lun,106)
      nl = 0
    END IF
    WRITE (lun,102) i, x(i), y(i)
    nl = nl + 1
  END DO
END IF

! Print the triangulation LTRI.

IF (nl > nlmax/2) THEN
  WRITE (lun,106)
  nl = 0
END IF
IF (nrow == 6) THEN
  WRITE (lun,103)
ELSE
  WRITE (lun,104)
END IF
nl = nl + 5
DO  k = 1,nt
  IF (nl >= nlmax) THEN
    WRITE (lun,106)
    nl = 0
  END IF
  WRITE (lun,105) k, ltri(1:nrow,k)
  nl = nl + 1
END DO

! Print NB, NA, and NT (boundary nodes, arcs, and triangles).

nb = 2*n - nt - 2
na = nt + n - 1
IF (nl > nlmax-6) WRITE (lun,106)
WRITE (lun,107) nb, na, nt

! Print NCC and LCT.

3 WRITE (lun,108) ncc
IF (ncc > 0) WRITE (lun,109) lct(1:ncc)
RETURN

! Print formats:

100 FORMAT (///t25, 'TRIPACK (TRLIST) Output')
101 FORMAT (//t17, 'Node       X(Node)          Y(Node)'//)
102 FORMAT (t17, i4, 2E17.6)
103 FORMAT (//' Triangle        Vertices            Neighbors'/  &
            '    KT       N1     N2     N3    KT1    KT2    KT3'/)
104 FORMAT (//' Triangle        Vertices            Neighbors', t61,'Arcs'/  &
    '    KT       N1     N2     N3    KT1    KT2    KT3    KA1    KA2    KA3'/)
105 FORMAT ('  ', i4, '  ', 6('   ', i4), 3('  ', i5))
106 FORMAT (///)
107 FORMAT (/' NB = ', i4, ' Boundary Nodes     NA = ', i5,  &
            ' Arcs     NT = ', i5, ' Triangles')
108 FORMAT (/' NCC =', i3, ' Constraint Curves')
109 FORMAT (t11, 14I5)
110 FORMAT (//t12, '*** Invalid Parameter:  N =', i5,  &
            ', NROW =', i5, ', NT =', i5, ' ***')
END SUBROUTINE trlprt



SUBROUTINE trmesh (n, x, y, list, lptr, lend, lnew, near, next, dist, ier)

INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN OUT)  :: x(:)
REAL, INTENT(IN OUT)  :: y(:)
INTEGER, INTENT(OUT)  :: list(:)
INTEGER, INTENT(OUT)  :: lptr(:)
INTEGER, INTENT(OUT)  :: lend(:)
INTEGER, INTENT(OUT)  :: lnew
INTEGER, INTENT(OUT)  :: near(:)
INTEGER, INTENT(OUT)  :: next(:)
REAL, INTENT(OUT)     :: dist(:)
INTEGER, INTENT(OUT)  :: ier

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   06/28/98

!   This subroutine creates a Delaunay triangulation of a set of N
! arbitrarily distributed points in the plane referred to as nodes.
! The Delaunay triangulation is defined as a set of triangles with the
! following five properties:

!  1)  The triangle vertices are nodes.
!  2)  No triangle contains a node other than its vertices.
!  3)  The interiors of the triangles are pairwise disjoint.
!  4)  The union of triangles is the convex hull of the set of nodes
!        (the smallest convex set which contains the nodes).
!  5)  The interior of the circumcircle of each triangle contains no node.

! The first four properties define a triangulation, and the last property
! results in a triangulation which is as close as possible to equiangular in
! a certain sense and which is uniquely defined unless four or more nodes lie
! on a common circle.  This property makes the triangulation well-suited
! for solving closest point problems and for triangle-based interpolation.

!   The triangulation can be generalized to a constrained Delaunay
! triangulation by a call to Subroutine ADDCST.
! This allows for user-specified boundaries defining a non-convex and/or
! multiply connected region.

!   The algorithm for constructing the triangulation has expected time
! complexity O(N*log(N)) for most nodal distributions.  Also, since the
! algorithm proceeds by adding nodes incrementally, the triangulation may be
! updated with the addition (or deletion) of a node very efficiently.
! The adjacency information representing the triangulation is stored as a
! linked list requiring approximately 13N storage locations.

!   The following is a list of the software package modules
! which a user may wish to call directly:

!  ADDCST - Generalizes the Delaunay triangulation to allow
!             for user-specified constraints.

!  ADDNOD - Updates the triangulation by appending or inserting a new node.

!  AREAP  - Computes the area bounded by a closed polygonal curve such as the
!             boundary of the triangulation or of a constraint region.

!  BNODES - Returns an array containing the indexes of the boundary nodes
!             in counterclockwise order.
!             Counts of boundary nodes, triangles, and arcs are also returned.

!  CIRCUM - Computes the area, circumcenter, circumradius, and, optionally,
!             the aspect ratio of a triangle defined by user-specified
!             vertices.

!  DELARC - Deletes a boundary arc from the triangulation.

!  DELNOD - Updates the triangulation with the deletion of a node.

!  EDGE   - Forces a pair of nodes to be connected by an arc
!             in the triangulation.

!  GETNP  - Determines the ordered sequence of L closest nodes to a given
!             node, along with the associated distances.  The distance
!             between nodes is taken to be the length of the shortest
!             connecting path which intersects no constraint region.

!  INTSEC - Determines whether or not an arbitrary pair of
!             line segments share a common point.

!  JRAND  - Generates a uniformly distributed pseudo-random integer.

!  LEFT   - Locates a point relative to a line.

!  NEARND - Returns the index of the nearest node to an
!             arbitrary point, along with its squared distance.

!  STORE  - Forces a value to be stored in main memory so
!             that the precision of floating point numbers
!             in memory locations rather than registers is computed.

!  TRLIST - Converts the triangulation data structure to a triangle
!             list more suitable for use in a finite element code.

!  TRLPRT - Prints the triangle list created by Subroutine TRLIST.

!  TRMESH - Creates a Delaunay triangulation of a set of nodes.

!  TRMSHR - Creates a Delaunay triangulation (more efficiently than TRMESH)
!             of a set of nodes lying at the vertices of a (possibly skewed)
!             rectangular grid.

!  TRPLOT - Creates a level-2 Encapsulated Postscript (EPS)
!             file containing a triangulation plot.

!  TRPRNT - Prints the triangulation data structure and,
!             optionally, the nodal coordinates.


! On input:

!       N = Number of nodes in the triangulation.  N >= 3.

!       X,Y = Arrays of length N containing the Cartesian coordinates of
!             the nodes.  (X(K),Y(K)) is referred to as node K, and K is
!             referred to as a nodal index.  The first three nodes must not
!             be collinear.

! The above parameters are not altered by this routine.

!       LIST,LPTR = Arrays of length at least 6N-12.

!       LEND = Array of length at least N.

!       NEAR,NEXT,DIST = Work space arrays of length at least N.
!                        The space is used to efficiently determine the
!                        nearest triangulation node to each unprocessed
!                        node for use by ADDNOD.

! On output:

!       LIST = Set of nodal indexes which, along with LPTR, LEND, and LNEW,
!              define the triangulation as a set of N adjacency lists
!              -- counterclockwise-ordered sequences of neighboring nodes such
!              that the first and last neighbors of a boundary node are
!              boundary nodes (the first neighbor of an interior node is
!              arbitrary).  In order to distinguish between interior and
!              boundary nodes, the last neighbor of each boundary node is
!              represented by the negative of its index.

!       LPTR = Set of pointers (LIST indexes) in one-to-one correspondence
!              with the elements of LIST.
!              LIST(LPTR(I)) indexes the node which follows LIST(I) in
!              cyclical counterclockwise order (the first neighbor follows
!              the last neighbor).

!       LEND = Set of pointers to adjacency lists.  LEND(K) points to the last
!              neighbor of node K for K = 1,...,N.
!              Thus, LIST(LEND(K)) < 0 if and only if K is a boundary node.

!       LNEW = Pointer to the first empty location in LIST and LPTR
!              (list length plus one).  LIST, LPTR, LEND, and LNEW are not
!              altered if IER < 0, and are incomplete if IER > 0.

!       NEAR,NEXT,DIST = Garbage.

!       IER = Error indicator:
!             IER =  0 if no errors were encountered.
!             IER = -1 if N < 3 on input.
!             IER = -2 if the first three nodes are collinear.
!             IER = -4 if an error flag was returned by a call to SWAP
!                      in ADDNOD.  This is an internal error and should be
!                      reported to the programmer.
!             IER =  L if nodes L and M coincide for some M > L.
!                      The linked list represents a triangulation of nodes
!                      1 to M-1 in this case.

! Modules required by TRMESH:  ADDNOD, BDYADD, INSERT, INTADD, JRAND, LEFT,
!                                LSTPTR, STORE, SWAP, SWPTST, TRFIND

! Intrinsic function called by TRMESH:  ABS

!***********************************************************

INTEGER :: i, i0, j, k, km1, lcc(1), lp, lpl, ncc, nexti, nn
REAL    :: d, d1, d2, d3, eps

! Local parameters:

! D =        Squared distance from node K to node I
! D1,D2,D3 = Squared distances from node K to nodes 1, 2, and 3, respectively
! EPS =      Half the machine precision
! I,J =      Nodal indexes
! I0 =       Index of the node preceding I in a sequence of
!              unprocessed nodes:  I = NEXT(I0)
! K =        Index of node to be added and DO-loop index: K > 3
! KM1 =      K-1
! LCC(1) =   Dummy array
! LP =       LIST index (pointer) of a neighbor of K
! LPL =      Pointer to the last neighbor of K
! NCC =      Number of constraint curves
! NEXTI =    NEXT(I)
! NN =       Local copy of N
! SWTOL =    Tolerance for function SWPTST

nn = n
IF (nn < 3) THEN
  ier = -1
  RETURN
END IF

! Compute a tolerance for function SWPTST:  SWTOL = 10*(machine precision)

eps = 1.
1 eps = eps/2.
swtol = store(eps + 1.)
IF (swtol > 1.) GO TO 1
swtol = eps*20.

! Store the first triangle in the linked list.

IF ( .NOT. left(x(1), y(1), x(2), y(2), x(3), y(3)) ) THEN
  
!   The initial triangle is (3,2,1) = (2,1,3) = (1,3,2).
  
  list(1) = 3
  lptr(1) = 2
  list(2) = -2
  lptr(2) = 1
  lend(1) = 2
  
  list(3) = 1
  lptr(3) = 4
  list(4) = -3
  lptr(4) = 3
  lend(2) = 4
  
  list(5) = 2
  lptr(5) = 6
  list(6) = -1
  lptr(6) = 5
  lend(3) = 6
  
ELSE IF ( .NOT. left(x(2), y(2), x(1), y(1), x(3), y(3)) ) THEN
  
!   The initial triangle is (1,2,3).
  
  list(1) = 2
  lptr(1) = 2
  list(2) = -3
  lptr(2) = 1
  lend(1) = 2
  
  list(3) = 3
  lptr(3) = 4
  list(4) = -1
  lptr(4) = 3
  lend(2) = 4
  
  list(5) = 1
  lptr(5) = 6
  list(6) = -2
  lptr(6) = 5
  lend(3) = 6
  
ELSE
  
!   The first three nodes are collinear.
  
  ier = -2
  RETURN
END IF

! Initialize LNEW and test for N = 3.

lnew = 7
IF (nn == 3) THEN
  ier = 0
  RETURN
END IF

! A nearest-node data structure (NEAR, NEXT, and DIST) is
!   used to obtain an expected-time (N*log(N)) incremental
!   algorithm by enabling constant search time for locating
!   each new node in the triangulation.

! For each unprocessed node K, NEAR(K) is the index of the
!   triangulation node closest to K (used as the starting
!   point for the search in Subroutine TRFIND) and DIST(K)
!   is an increasing function of the distance between nodes
!   K and NEAR(K).

! Since it is necessary to efficiently find the subset of
!   unprocessed nodes associated with each triangulation
!   node J (those that have J as their NEAR entries), the
!   subsets are stored in NEAR and NEXT as follows:  for
!   each node J in the triangulation, I = NEAR(J) is the
!   first unprocessed node in J's set (with I = 0 if the
!   set is empty), L = NEXT(I) (if I > 0) is the second,
!   NEXT(L) (if L > 0) is the third, etc.  The nodes in each
!   set are initially ordered by increasing indexes (which
!   maximizes efficiency) but that ordering is not main-
!   tained as the data structure is updated.

! Initialize the data structure for the single triangle.

near(1) = 0
near(2) = 0
near(3) = 0
DO  k = nn,4,-1
  d1 = (x(k)-x(1))**2 + (y(k)-y(1))**2
  d2 = (x(k)-x(2))**2 + (y(k)-y(2))**2
  d3 = (x(k)-x(3))**2 + (y(k)-y(3))**2
  IF (d1 <= d2 .AND. d1 <= d3) THEN
    near(k) = 1
    dist(k) = d1
    next(k) = near(1)
    near(1) = k
  ELSE IF (d2 <= d1 .AND. d2 <= d3) THEN
    near(k) = 2
    dist(k) = d2
    next(k) = near(2)
    near(2) = k
  ELSE
    near(k) = 3
    dist(k) = d3
    next(k) = near(3)
    near(3) = k
  END IF
END DO

! Add the remaining nodes.  Parameters for ADDNOD are as follows:

!   K = Index of the node to be added.
!   NEAR(K) = Index of the starting node for the search in TRFIND.
!   NCC = Number of constraint curves.
!   LCC = Dummy array (since NCC = 0).
!   KM1 = Number of nodes in the triangulation.

ncc = 0
DO  k = 4,nn
  km1 = k-1
  CALL addnod (k, x(k), y(k), near(k), ncc, lcc, km1, x, y, list, lptr,  &
               lend, lnew, ier)
  IF (ier /= 0) RETURN
  
! Remove K from the set of unprocessed nodes associated with NEAR(K).
  
  i = near(k)
  IF (near(i) == k) THEN
    near(i) = next(k)
  ELSE
    i = near(i)
    3 i0 = i
    i = next(i0)
    IF (i /= k) GO TO 3
    next(i0) = next(k)
  END IF
  near(k) = 0
  
! Loop on neighbors J of node K.
  
  lpl = lend(k)
  lp = lpl
  4   lp = lptr(lp)
  j = ABS(list(lp))
  
! Loop on elements I in the sequence of unprocessed nodes associated with J:
!   K is a candidate for replacing J as the nearest triangulation node to I.
!   The next value of I in the sequence, NEXT(I), must be saved before I
!   is moved because it is altered by adding I to K's set.
  
  i = near(j)
  5 IF (i == 0) GO TO 6
  nexti = next(i)
  
! Test for the distance from I to K less than the distance
!   from I to J.
  
  d = (x(k)-x(i))**2 + (y(k)-y(i))**2
  IF (d < dist(i)) THEN
    
! Replace J by K as the nearest triangulation node to I:
!   update NEAR(I) and DIST(I), and remove I from J's set
!   of unprocessed nodes and add it to K's set.
    
    near(i) = k
    dist(i) = d
    IF (i == near(j)) THEN
      near(j) = nexti
    ELSE
      next(i0) = nexti
    END IF
    next(i) = near(k)
    near(k) = i
  ELSE
    i0 = i
  END IF
  
! Bottom of loop on I.
  
  i = nexti
  GO TO 5
  
! Bottom of loop on neighbors J.
  
  6 IF (lp /= lpl) GO TO 4
END DO
RETURN
END SUBROUTINE trmesh



SUBROUTINE trmshr (n, nx, x, y, nit, list, lptr, lend, lnew, ier)

INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN)      :: nx
REAL, INTENT(IN)         :: x(:)
REAL, INTENT(IN)         :: y(:)
INTEGER, INTENT(IN OUT)  :: nit
INTEGER, INTENT(OUT)     :: list(:)
INTEGER, INTENT(OUT)     :: lptr(:)
INTEGER, INTENT(OUT)     :: lend(:)
INTEGER, INTENT(OUT)     :: lnew
INTEGER, INTENT(OUT)     :: ier

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   06/27/98

!   This subroutine creates a Delaunay triangulation of a set of N nodes in
! the plane, where the nodes are the vertices of an NX by NY skewed
! rectangular grid with the natural ordering.   Thus, N = NX*NY, and the nodes
! are ordered from left to right beginning at the top row so that adjacent
! nodes have indexes which differ by 1 in the x-direction and by NX in the
! y-direction.  A skewed rectangular grid is defined as one in which each grid
! cell is a strictly convex quadrilateral (and is thus the convex hull of its
! four vertices).  Equivalently, any transformation from a rectangle to a grid
! cell which is bilinear in both components has an invertible Jacobian.

!   If the nodes are not distributed and ordered as defined above,
! Subroutine TRMESH must be called in place of this routine.
! Refer to Subroutine ADDCST for the treatment of constraints.

!   The first phase of the algorithm consists of constructing a triangulation
! by choosing a diagonal arc in each grid cell.  If NIT = 0, all diagonals
! connect lower left to upper right corners and no error checking or
! additional computation is performed.  Otherwise, each diagonal arc is
! chosen to be locally optimal, and boundary arcs are added where necessary
! in order to cover the convex hull of the nodes.  (This is the first
! iteration.)  If NIT > 1 and no error was detected, the triangulation is then
! optimized by a sequence of up to NIT-1 iterations in which interior arcs of
! the triangulation are tested and swapped if appropriate.  The algorithm
! terminates when an iteration results in no swaps and/or when the allowable
! number of iterations has been performed.  NIT = 0 is sufficient to produce
! a Delaunay triangulation if the original grid is actually rectangular, and
! NIT = 1 is sufficient if it is close to rectangular.  Note, however, that
! the ordering and distribution of nodes is not checked for validity in the
! case NIT = 0, and the triangulation will not be valid unless the rectangular
! grid covers the convex hull of the nodes.

! On input:

!       N = Number of nodes in the grid.  N = NX*NY for some NY >= 2.

!       NX = Number of grid points in the x-direction.  NX >= 2.

!       X,Y = Arrays of length N containing coordinates of the nodes with the
!             ordering and distribution defined in the header comments above.
!             (X(K),Y(K)) is referred to as node K.

! The above parameters are not altered by this routine.

!       NIT = Nonnegative integer specifying the maximum number of iterations
!             to be employed.  Refer to the header comments above.

!       LIST,LPTR = Arrays of length at least 6N-12.

!       LEND = Array of length at least N.

! On output:

!       NIT = Number of iterations employed.

!       LIST,LPTR,LEND,LNEW = Data structure defining the
!                             triangulation.  Refer to Subroutine TRMESH.

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = K if the grid element with upper left
!                     corner at node K is not a strictly
!                     convex quadrilateral.  The algorithm
!                     is terminated when the first such
!                     occurrence is detected.  Note that
!                     this test is not performed if NIT = 0
!                     on input.
!             IER = -1 if N, NX, or NIT is outside its valid range on input.
!             IER = -2 if NIT > 1 on input, and the optimi-
!                      zation loop failed to converge within
!                      the allowable number of iterations.
!                      The triangulation is valid but not
!                      optimal in this case.

! Modules required by TRMSHR:  INSERT, LEFT, LSTPTR, NBCNT,
!                                STORE, SWAP, SWPTST

! Intrinsic function called by TRMSHR:  ABS

!***********************************************************

INTEGER :: i, iter, j, k, kp1, lp, lpf, lpk, lpl, lpp,  &
           m1, m2, m3, m4, maxit, n0, n1, n2, n3, n4, ni, nj, nm1, nn, nnb
LOGICAL :: tst
REAL    :: eps

! Store local variables and test for errors in input parameters.

ni = nx
nj = n/ni
nn = ni*nj
maxit = nit
nit = 0
IF (n /= nn .OR. nj < 2 .OR. ni < 2  .OR. maxit < 0) THEN
  ier = -1
  RETURN
END IF
ier = 0

! Compute a tolerance for function SWPTST:  SWTOL = 10*(machine precision)

eps = 1.
1 eps = eps/2.
swtol = store(eps + 1.)
IF (swtol > 1.) GO TO 1
swtol = eps*20.

! Loop on grid points (I,J) corresponding to nodes K = (J-1)*NI + I.
!   TST = TRUE iff diagonals are to be chosen by the swap test.  M1, M2, M3,
!   and M4 are the slopes (-1, 0, or 1) of the diagonals in quadrants 1 to 4
!   (counterclockwise beginning with the upper right)
!   for a coordinate system with origin at node K.

tst = maxit > 0
m1 = 0
m4 = 0
lp = 0
kp1 = 1
DO  j = 1,nj
  DO  i = 1,ni
    m2 = m1
    m3 = m4
    k = kp1
    kp1 = k + 1
    lpf = lp + 1
    IF (j == nj .AND. i /= ni) GO TO 2
    IF (i /= 1) THEN
      IF (j /= 1) THEN
        
!   K is not in the top row, leftmost column, or bottom row (unless K is the
!     lower right corner).  Take the first neighbor to be the node above K.
        
        lp = lp + 1
        list(lp) = k - ni
        lptr(lp) = lp + 1
        IF (m2 <= 0) THEN
          lp = lp + 1
          list(lp) = k - 1 - ni
          lptr(lp) = lp + 1
        END IF
      END IF
      
!   K is not in the leftmost column.  The next (or first)
!     neighbor is to the left of K.
      
      lp = lp + 1
      list(lp) = k - 1
      lptr(lp) = lp + 1
      IF (j == nj) GO TO 3
      IF (m3 >= 0) THEN
        lp = lp + 1
        list(lp) = k - 1 + ni
        lptr(lp) = lp + 1
      END IF
    END IF
    
!   K is not in the bottom row.  The next (or first) neighbor is below K.
    
    lp = lp + 1
    list(lp) = k + ni
    lptr(lp) = lp + 1
    
!   Test for a negative diagonal in quadrant 4 unless K is in the rightmost
!     column.  The quadrilateral associated with the quadrant is tested for
!     strict convexity unless NIT = 0 on input.
    
    IF (i == ni) GO TO 3
    m4 = 1
    IF (.NOT. tst) GO TO 2
    IF ( left(x(kp1),y(kp1),x(k+ni),y(k+ni),x(k),y(k))  &
        .OR. left(x(k),y(k),x(kp1+ni),y(kp1+ni), x(k+ni),y(k+ni))  &
        .OR. left(x(k+ni),y(k+ni),x(kp1),y(kp1), x(kp1+ni),y(kp1+ni))  &
        .OR. left(x(kp1+ni),y(kp1+ni),x(k),y(k),  &
        x(kp1),y(kp1)) ) GO TO 12
    IF ( swptst(kp1,k+ni,k,kp1+ni,x,y) ) GO TO 2
    m4 = -1
    lp = lp + 1
    list(lp) = kp1 + ni
    lptr(lp) = lp + 1
    
!   The next (or first) neighbor is to the right of K.
    
    2 lp = lp + 1
    list(lp) = kp1
    lptr(lp) = lp + 1
    
!   Test for a positive diagonal in quadrant 1 (the neighbor of K-NI
!     which follows K is not K+1) unless K is in the top row.
    
    IF (j == 1) GO TO 3
    IF (tst) THEN
      m1 = -1
      lpk = lstptr(lend(k-ni),k,list,lptr)
      lpk = lptr(lpk)
      IF (list(lpk) /= kp1) THEN
        m1 = 1
        lp = lp + 1
        list(lp) = kp1 - ni
        lptr(lp) = lp + 1
      END IF
    END IF
    
!   If K is in the leftmost column (and not the top row) or in the bottom row
!     (and not the rightmost column), then the next neighbor is the node
!     above K.
    
    IF (i /= 1 .AND. j /= nj) GO TO 4
    lp = lp + 1
    list(lp) = k - ni
    lptr(lp) = lp + 1
    IF (i == 1) GO TO 3
    
!   K is on the bottom row (and not the leftmost or rightmost column).
    
    IF (m2 <= 0) THEN
      lp = lp + 1
      list(lp) = k - 1 - ni
      lptr(lp) = lp + 1
    END IF
    lp = lp + 1
    list(lp) = k - 1
    lptr(lp) = lp + 1
    
!   K is a boundary node.
    
    3 list(lp) = -list(lp)
    
!   Bottom of loop.  Store LEND and correct LPTR(LP).
!     LPF and LP point to the first and last neighbors of K.
    
    4 lend(k) = lp
    lptr(lp) = lpf
  END DO
END DO

! Store LNEW, and terminate the algorithm if NIT = 0 on input.

lnew = lp + 1
IF (maxit == 0) RETURN

! Add boundary arcs where necessary in order to cover the convex hull of the
!   nodes.  N1, N2, and N3 are consecutive boundary nodes in counterclockwise
!   order, and N0 is the starting point for each loop around the boundary.

n0 = 1
n1 = n0
n2 = ni + 1

!   TST is set to TRUE if an arc is added.  The boundary
!     loop is repeated until a traversal results in no added arcs.

7 tst = .false.

!   Top of boundary loop.  Set N3 to the first neighbor of
!     N2, and test for N3 LEFT N1 -> N2.

8   lpl = lend(n2)
lp = lptr(lpl)
n3 = list(lp)
IF ( left(x(n1),y(n1),x(n2),y(n2),x(n3),y(n3)) ) n1 = n2
IF (n1 /= n2) THEN
  
!   Add the boundary arc N1-N3.  If N0 = N2, the starting
!     point is changed to N3, since N2 will be removed from
!     the boundary.  N3 is inserted as the first neighbor of
!     N1, N2 is changed to an interior node, and N1 is
!     inserted as the last neighbor of N3.
  
  tst = .true.
  IF (n2 == n0) n0 = n3
  lp = lend(n1)
  CALL insert (n3,lp, list,lptr,lnew )
  list(lpl) = -list(lpl)
  lp = lend(n3)
  list(lp) = n2
  CALL insert (-n1,lp, list,lptr,lnew )
  lend(n3) = lnew - 1
END IF

!   Bottom of loops.  Test for termination.

n2 = n3
IF (n1 /= n0) GO TO 8
IF (tst) GO TO 7

! Terminate the algorithm if NIT = 1 on input.

nit = 1
IF (maxit == 1) RETURN

! Optimize the triangulation by applying the swap test and appropriate swaps
!   to the interior arcs.  The loop is repeated until no swaps are performed
!   or MAXIT iterations have been applied.  ITER is the current iteration,
!   and TST is set to TRUE if a swap occurs.

iter = 1
nm1 = nn - 1
9 iter = iter + 1
tst = .false.

!   Loop on interior arcs N1-N2, where N2 > N1 and
!     (N1,N2,N3) and (N2,N1,N4) are adjacent triangles.

!   Top of loop on nodes N1.

loop11:  DO  n1 = 1,nm1
  lpl = lend(n1)
  n4 = list(lpl)
  lpf = lptr(lpl)
  n2 = list(lpf)
  lp = lptr(lpf)
  n3 = list(lp)
  nnb = nbcnt(lpl,lptr)
  
!   Top of loop on neighbors N2 of N1.  NNB is the number of neighbors of N1.
  
  DO  i = 1,nnb
    
!   Bypass the swap test if N1 is a boundary node and N2 is
!     the first neighbor (N4 < 0), N2 < N1, or N1-N2 is a
!     diagonal arc (already locally optimal) when ITER = 2.
    
    IF ( n4 > 0 .AND. n2 > n1 .AND. &
        (iter /= 2 .OR. ABS(n1+ni-n2) /= 1) ) THEN
    IF (swptst(n3, n4, n1, n2, x, y) ) THEN
      
!   Swap diagonal N1-N2 for N3-N4, set TST to TRUE, and set
!     N2 to N4 (the neighbor preceding N3).
      
      CALL swap (n3, n4, n1, n2, list, lptr, lend, lpp)
      IF (lpp /= 0) THEN
        tst = .true.
        n2 = n4
      END IF
    END IF
  END IF
  
!   Bottom of neighbor loop.
  
  IF (list(lpl) == -n3) CYCLE loop11
  n4 = n2
  n2 = n3
  lp = lstptr(lpl,n2,list,lptr)
  lp = lptr(lp)
  n3 = ABS(list(lp))
END DO
END DO loop11

!   Test for termination.

IF (tst .AND. iter < maxit) GO TO 9
nit = iter
IF (tst) ier = -2
RETURN

! Invalid grid cell encountered.

12 ier = k
RETURN
END SUBROUTINE trmshr



SUBROUTINE trplot (lun,pltsiz,wx1,wx2,wy1,wy2,ncc,lcc,  &
                   n,x,y,list,lptr,lend,title, numbr, ier)

INTEGER, INTENT(IN)            :: lun
REAL, INTENT(IN)               :: pltsiz
REAL, INTENT(IN)               :: wx1
REAL, INTENT(IN)               :: wx2
REAL, INTENT(IN)               :: wy1
REAL, INTENT(IN)               :: wy2
INTEGER, INTENT(IN)            :: ncc
INTEGER, INTENT(IN)            :: lcc(:)
INTEGER, INTENT(IN)            :: n
REAL, INTENT(IN)               :: x(n)
REAL, INTENT(IN)               :: y(n)
INTEGER, INTENT(IN OUT)        :: list(:)
INTEGER, INTENT(IN)            :: lptr(:)
INTEGER, INTENT(IN)            :: lend(:)
CHARACTER (LEN=*), INTENT(IN)  :: title
LOGICAL, INTENT(IN)            :: numbr
INTEGER, INTENT(OUT)           :: ier

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/15/98

!   This subroutine creates a level-2 Encapsulated Postscript (EPS) file
! containing a triangulation plot.


! On input:

!       LUN = Logical unit number in the range 0 to 99.
!             The unit should be opened with an appropriate file name before
!             the call to this routine.

!       PLTSIZ = Plot size in inches.  The window is mapped, with aspect ratio
!                preserved, to a rectangular viewport with maximum side-length
!                equal to .88*PLTSIZ (leaving room for labels outside the
!                viewport).  The viewport is centered on the 8.5 by 11 inch
!                page, and its boundary is drawn.  1.0 <= PLTSIZ <= 8.5.

!       WX1,WX2,WY1,WY2 = Parameters defining a rectangular window against
!                         which the triangulation is clipped.  (Only the
!                         portion of the triangulation that lies in the window
!                         is drawn.)
!                         (WX1,WY1) and (WX2,WY2) are the lower left and upper
!                         right corners, respectively.  WX1 < WX2 and
!                         WY1 < WY2.

!       NCC = Number of constraint curves.  Refer to Subroutine ADDCST.
!             NCC >= 0.

!       LCC = Array of length NCC (or dummy parameter if
!             NCC = 0) containing the index of the first
!             node of constraint I in LCC(I).  For I = 1 to
!             NCC, LCC(I+1)-LCC(I) >= 3, where LCC(NCC+1) = N+1.

!       N = Number of nodes in the triangulation.  N >= 3.

!       X,Y = Arrays of length N containing the coordinates of the nodes with
!             non-constraint nodes in the first LCC(1)-1 locations.

!       LIST,LPTR,LEND = Data structure defining the triangulation.
!                        Refer to Subroutine TRMESH.

!       TITLE = Type CHARACTER variable or constant containing a string to be
!               centered above the plot.
!               The string must be enclosed in parentheses; i.e., the first
!               and last characters must be '(' and ')', respectively, but
!               these are not displayed.  TITLE may have at most 80 characters
!               including the parentheses.

!       NUMBR = Option indicator:  If NUMBR = TRUE, the
!               nodal indexes are plotted next to the nodes.

! Input parameters are not altered by this routine.

! On output:

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if LUN, PLTSIZ, NCC, or N is outside its valid range.
!                     LCC is not tested for validity.
!             IER = 2 if WX1 >= WX2 or WY1 >= WY2.
!             IER = 3 if an error was encountered in writing to unit LUN.

!   Various plotting options can be controlled by altering
! the data statement below.

! Modules required by TRPLOT:  None

! Intrinsic functions called by TRPLOT:  ABS, CHAR, NINT, REAL

!***********************************************************

INTEGER  :: i, ifrst, ih, ilast, ipx1, ipx2, ipy1, ipy2,  &
            iw, lp, lpl, n0, n0bak, n0for, n1, nls
LOGICAL  :: cnstr, pass1
REAL     :: dx, dy, r, sfx, sfy, t, tx, ty, x0, y0
LOGICAL, PARAMETER  :: annot = .TRUE.
REAL, PARAMETER     :: dashl = 4.0, fsizn = 10.0, fsizt = 16.0

! Local parameters:

! ANNOT =     Logical variable with value TRUE iff the plot
!               is to be annotated with the values of WX1, WX2, WY1, and WY2
! CNSTR       Logical variable used to flag constraint arcs:
!               TRUE iff N0-N1 lies in a constraint region
! DASHL =     Length (in points, at 72 points per inch) of dashes and spaces
!               in a dashed line pattern used for drawing constraint arcs
! DX =        Window width WX2-WX1
! DY =        Window height WY2-WY1
! FSIZN =     Font size in points for labeling nodes with
!               their indexes if NUMBR = TRUE
! FSIZT =     Font size in points for the title (and
!               annotation if ANNOT = TRUE)
! I =         Constraint index (1 to NCC)
! IFRST =     Index of the first node in constraint I
! IH =        Height of the viewport in points
! ILAST =     Index of the last node in constraint I
! IPX1,IPY1 = X and y coordinates (in points) of the lower
!               left corner of the bounding box or viewport
! IPX2,IPY2 = X and y coordinates (in points) of the upper
!               right corner of the bounding box or viewport
! IW =        Width of the viewport in points
! LP =        LIST index (pointer)
! LPL =       Pointer to the last neighbor of N0
! N0 =        Nodal index and DO-loop index
! N0BAK =     Predecessor of N0 in a constraint curve
!               (sequence of adjacent constraint nodes)
! N0FOR =     Successor to N0 in a constraint curve
! N1 =        Index of a neighbor of N0
! NLS =       Index of the last non-constraint node
! PASS1 =     Logical variable used to flag the first pass
!               through the constraint nodes
! R =         Aspect ratio DX/DY
! SFX,SFY =   Scale factors for mapping world coordinates
!               (window coordinates in [WX1,WX2] X [WY1,WY2])
!               to viewport coordinates in [IPX1,IPX2] X [IPY1,IPY2]
! T =         Temporary variable
! TX,TY =     Translation vector for mapping world coordinates to
!               viewport coordinates
! X0,Y0 =     X(N0),Y(N0) or label location


! Test for error 1, and set NLS to the last non-constraint node.

IF (lun < 0 .OR. lun > 99  .OR. pltsiz < 1.0 .OR. pltsiz > 8.5 .OR. &
    ncc < 0 .OR. n < 3) GO TO 11
nls = n
IF (ncc > 0) nls = lcc(1)-1

! Compute the aspect ratio of the window.

dx = wx2 - wx1
dy = wy2 - wy1
IF (dx <= 0.0 .OR. dy <= 0.0) GO TO 12
r = dx/dy

! Compute the lower left (IPX1,IPY1) and upper right
!   (IPX2,IPY2) corner coordinates of the bounding box.
!   The coordinates, specified in default user space units
!   (points, at 72 points/inch with origin at the lower
!   left corner of the page), are chosen to preserve the
!   aspect ratio R, and to center the plot on the 8.5 by 11
!   inch page.  The center of the page is (306,396), and
!   T = PLTSIZ/2 in points.

t = 36.0*pltsiz
IF (r >= 1.0) THEN
  ipx1 = 306 - nint(t)
  ipx2 = 306 + nint(t)
  ipy1 = 396 - nint(t/r)
  ipy2 = 396 + nint(t/r)
ELSE
  ipx1 = 306 - nint(t*r)
  ipx2 = 306 + nint(t*r)
  ipy1 = 396 - nint(t)
  ipy2 = 396 + nint(t)
END IF

! Output header comments.

WRITE (lun,100,ERR=13) ipx1, ipy1, ipx2, ipy2
100 FORMAT ('%!PS-Adobe-3.0 EPSF-3.0'/ '%%BoundingBox:',4I4/  &
    '%%Title:  Triangulation'/ '%%Creator:  TRIPACK'/  &
    '%%EndComments')

! Set (IPX1,IPY1) and (IPX2,IPY2) to the corner coordinates
!   of a viewport obtained by shrinking the bounding box by
!   12% in each dimension.

iw = nint(0.88*REAL(ipx2-ipx1))
ih = nint(0.88*REAL(ipy2-ipy1))
ipx1 = 306 - iw/2
ipx2 = 306 + iw/2
ipy1 = 396 - ih/2
ipy2 = 396 + ih/2

! Set the line thickness to 2 points, and draw the viewport boundary.

t = 2.0
WRITE (lun,110,ERR=13) t
WRITE (lun,120,ERR=13) ipx1, ipy1
WRITE (lun,130,ERR=13) ipx1, ipy2
WRITE (lun,130,ERR=13) ipx2, ipy2
WRITE (lun,130,ERR=13) ipx2, ipy1
WRITE (lun,140,ERR=13)
WRITE (lun,150,ERR=13)
110 FORMAT (f12.6,' setlinewidth')
120 FORMAT (2I4,' moveto')
130 FORMAT (2I4,' lineto')
140 FORMAT ('closepath')
150 FORMAT ('stroke')

! Set up a mapping from the window to the viewport.

sfx = REAL(iw)/dx
sfy = REAL(ih)/dy
tx = ipx1 - sfx*wx1
ty = ipy1 - sfy*wy1
WRITE (lun,160,ERR=13) tx, ty, sfx, sfy
160 FORMAT (2F12.6, ' translate'/ 2F12.6, ' scale')

! The line thickness (believe it or fucking not) must be
!   changed to reflect the new scaling which is applied to
!   all subsequent output.  Set it to 1.0 point.

t = 2.0/(sfx+sfy)
WRITE (lun,110,ERR=13) t

! Save the current graphics state, and set the clip path to
!   the boundary of the window.

WRITE (lun,170,ERR=13)
WRITE (lun,180,ERR=13) wx1, wy1
WRITE (lun,190,ERR=13) wx2, wy1
WRITE (lun,190,ERR=13) wx2, wy2
WRITE (lun,190,ERR=13) wx1, wy2
WRITE (lun,200,ERR=13)
170 FORMAT ('gsave')
180 FORMAT (2F12.6,' moveto')
190 FORMAT (2F12.6,' lineto')
200 FORMAT ('closepath clip newpath')

! Draw the edges N0->N1, where N1 > N0, beginning with a loop on
!   non-constraint nodes N0.  LPL points to the last neighbor of N0.

DO  n0 = 1,nls
  x0 = x(n0)
  y0 = y(n0)
  lpl = lend(n0)
  lp = lpl
  
!   Loop on neighbors N1 of N0.
  
  2 lp = lptr(lp)
  n1 = ABS(list(lp))
  IF (n1 > n0) THEN
    
!   Add the edge to the path.
    
    WRITE (lun,210,ERR=13) x0, y0, x(n1), y(n1)
    210 FORMAT (2F12.6,' moveto',2F12.6,' lineto')
  END IF
  IF (lp /= lpl) GO TO 2
END DO

! Loop through the constraint nodes twice.  The non-
!   constraint arcs incident on constraint nodes are
!   drawn (with solid lines) on the first pass, and the
!   constraint arcs (both boundary and interior, if any)
!   are drawn (with dashed lines) on the second pass.

pass1 = .true.

! Loop on constraint nodes N0 with (N0BAK,N0,N0FOR) a sub-
!   sequence of constraint I.  The outer loop is on
!   constraints I with first and last nodes IFRST and ILAST.

4 ifrst = n+1
DO  i = ncc,1,-1
  ilast = ifrst - 1
  ifrst = lcc(i)
  n0bak = ilast
  DO  n0 = ifrst,ilast
    n0for = n0 + 1
    IF (n0 == ilast) n0for = ifrst
    lpl = lend(n0)
    x0 = x(n0)
    y0 = y(n0)
    lp = lpl
    
!   Loop on neighbors N1 of N0.  CNSTR = TRUE iff N0-N1 is a
!     constraint arc.
    
!   Initialize CNSTR to TRUE iff the first neighbor of N0
!     strictly follows N0FOR and precedes or coincides with
!     N0BAK (in counterclockwise order).
    
    5 lp = lptr(lp)
    n1 = ABS(list(lp))
    IF (n1 /= n0for .AND. n1 /= n0bak) GO TO 5
    cnstr = n1 == n0bak
    lp = lpl
    
!   Loop on neighbors N1 of N0.  Update CNSTR and test for
!     N1 > N0.
    
    6 lp = lptr(lp)
    n1 = ABS(list(lp))
    IF (n1 == n0for) cnstr = .true.
    IF (n1 > n0) THEN
      
!   Draw the edge iff (PASS1=TRUE and CNSTR=FALSE) or
!     (PASS1=FALSE and CNSTR=TRUE); i.e., CNSTR and PASS1
!     have opposite values.
      
      IF (cnstr .NEQV. pass1) WRITE (lun,210,ERR=13) x0, y0, x(n1), y(n1)
    END IF
    IF (n1 == n0bak) cnstr = .false.
    
!   Bottom of loops.
    
    IF (lp /= lpl) GO TO 6
    n0bak = n0
  END DO
END DO
IF (pass1) THEN
  
! End of first pass:  paint the path and change to dashed
!   lines for subsequent drawing.  Since the scale factors
!   are applied to everything, the dash length must be
!   specified in world coordinates.
  
  pass1 = .false.
  WRITE (lun,150,ERR=13)
  t = dashl*2.0/(sfx+sfy)
  WRITE (lun,220,ERR=13) t
  220 FORMAT ('[',f12.6,'] 0 setdash')
  GO TO 4
END IF

! Paint the path and restore the saved graphics state (with
!   no clip path).

WRITE (lun,150,ERR=13)
WRITE (lun,230,ERR=13)
230 FORMAT ('grestore')
IF (numbr) THEN
  
! Nodes in the window are to be labeled with their indexes.
!   Convert FSIZN from points to world coordinates, and
!   output the commands to select a font and scale it.
  
  t = fsizn*2.0/(sfx+sfy)
  WRITE (lun,240,ERR=13) t
  240 FORMAT ('/Helvetica findfont'/ f12.6,' scalefont setfont')
  
!   Loop on nodes N0 with coordinates (X0,Y0).
  
  DO  n0 = 1,n
    x0 = x(n0)
    y0 = y(n0)
    IF (x0 < wx1 .OR. x0 > wx2  .OR. y0 < wy1 .OR. y0 > wy2) CYCLE
    
!   Move to (X0,Y0), and draw the label N0.  The first char-
!     acter will have its lower left corner about one
!     character width to the right of the nodal position.
    
    WRITE (lun,180,ERR=13) x0, y0
    WRITE (lun,250,ERR=13) n0
    250 FORMAT ('(',i3,') show')
  END DO
END IF

! Convert FSIZT from points to world coordinates, and output
!   the commands to select a font and scale it.

t = fsizt*2.0/(sfx+sfy)
WRITE (lun,240,ERR=13) t

! Display TITLE centered above the plot:

y0 = wy2 + 3.0*t
WRITE (lun,260,ERR=13) title, (wx1+wx2)/2.0, y0
260 FORMAT (a80/'  stringwidth pop 2 div neg ',f12.6, ' add ',f12.6,' moveto')
WRITE (lun,270,ERR=13) title
270 FORMAT (a80/'  show')

IF (annot) THEN
  
! Display the window extrema below the plot.
  
  x0 = wx1
  y0 = wy1 - 100.0/(sfx+sfy)
  WRITE (lun,180,ERR=13) x0, y0
  WRITE (lun,280,ERR=13) wx1, wx2
  y0 = y0 - 2.0*t
  WRITE (lun,290,ERR=13) x0, y0, wy1, wy2
  280 FORMAT ('(Window:   WX1 = ',e9.3,',   WX2 = ',e9.3, ') show')
  290 FORMAT ('(Window:  ) stringwidth pop ',f12.6,' add',  &
              f12.6,' moveto'/ '( WY1 = ',e9.3,',   WY2 = ',e9.3,') show')
END IF

! Paint the path and output the showpage command and
!   end-of-file indicator.

WRITE (lun,300,ERR=13)
300 FORMAT ('stroke'/ 'showpage'/ '%%EOF')

! HP's interpreters require a one-byte End-of-PostScript-Job
!   indicator (to eliminate a timeout error message):
!   ASCII 4.

WRITE (lun,310,ERR=13) CHAR(4)
310 FORMAT (a1)

! No error encountered.

ier = 0
RETURN

! Invalid input parameter.

11 ier = 1
RETURN

! DX or DY is not positive.

12 ier = 2
RETURN

! Error writing to unit LUN.

13 ier = 3
RETURN
END SUBROUTINE trplot



SUBROUTINE trprnt (ncc, lcc, n, x, y, list, lptr, lend, lout, prntx)

INTEGER, INTENT(IN)      :: ncc
INTEGER, INTENT(IN OUT)  :: lcc(:)
INTEGER, INTENT(IN)      :: n
REAL, INTENT(IN OUT)     :: x(n)
REAL, INTENT(IN OUT)     :: y(n)
INTEGER, INTENT(IN)      :: list(:)
INTEGER, INTENT(IN)      :: lptr(:)
INTEGER, INTENT(IN)      :: lend(:)
INTEGER, INTENT(IN)      :: lout
LOGICAL, INTENT(IN)      :: prntx

!***********************************************************

!                                               From TRIPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/30/98

!   Given a triangulation of a set of points in the plane,
! this subroutine prints the adjacency lists and, option-
! ally, the nodal coordinates on logical unit LOUT.  The
! list of neighbors of a boundary node is followed by index
! 0.  The numbers of boundary nodes, triangles, and arcs,
! and the constraint curve starting indexes, if any, are
! also printed.


! On input:

!       NCC = Number of constraints.

!       LCC = List of constraint curve starting indexes (or
!             dummy array of length 1 if NCC = 0).

!       N = Number of nodes in the triangulation.
!           3 <= N <= 9999.

!       X,Y = Arrays of length N containing the coordinates
!             of the nodes in the triangulation -- not used
!             unless PRNTX = TRUE.

!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to Subroutine
!                        TRMESH.

!       LOUT = Logical unit number for output.  0 <= LOUT
!              <= 99.  Output is printed on unit 6 if LOUT
!              is outside its valid range on input.

!       PRNTX = Logical variable with value TRUE if and only
!               if X and Y are to be printed (to 6 decimal places).

! None of the parameters are altered by this routine.

! Modules required by TRPRNT:  None

!***********************************************************

INTEGER  :: inc, k, lp, lpl, lun, na, nabor(100), nb, nd, nl, node, nn, nt
INTEGER, PARAMETER  :: nmax = 9999, nlmax = 60

nn = n
lun = lout
IF (lun < 0 .OR. lun > 99) lun = 6

! Print a heading and test the range of N.

WRITE (lun,100) nn
IF (nn < 3 .OR. nn > nmax) THEN
  
! N is outside its valid range.
  
  WRITE (lun,110)
  GO TO 5
END IF

! Initialize NL (the number of lines printed on the current
!   page) and NB (the number of boundary nodes encountered).

nl = 6
nb = 0
IF (.NOT. prntx) THEN
  
! Print LIST only.  K is the number of neighbors of NODE
!   which are stored in NABOR.
  
  WRITE (lun,101)
  DO  node = 1,nn
    lpl = lend(node)
    lp = lpl
    k = 0
    
    1 k = k + 1
    lp = lptr(lp)
    nd = list(lp)
    nabor(k) = nd
    IF (lp /= lpl) GO TO 1
    IF (nd <= 0) THEN
      
!   NODE is a boundary node.  Correct the sign of the last
!     neighbor, add 0 to the end of the list, and increment NB.
      
      nabor(k) = -nd
      k = k + 1
      nabor(k) = 0
      nb = nb + 1
    END IF
    
!   Increment NL and print the list of neighbors.
    
    inc = (k-1)/14 + 2
    nl = nl + inc
    IF (nl > nlmax) THEN
      WRITE (lun,106)
      nl = inc
    END IF
    WRITE (lun,103) node, nabor(1:k)
    IF (k /= 14) WRITE (lun,105)
  END DO
ELSE
  
! Print X, Y, and LIST.
  
  WRITE (lun,102)
  DO  node = 1,nn
    lpl = lend(node)
    lp = lpl
    k = 0
    3     k = k + 1
    lp = lptr(lp)
    nd = list(lp)
    nabor(k) = nd
    IF (lp /= lpl) GO TO 3
    IF (nd <= 0) THEN
      
!   NODE is a boundary node.
      
      nabor(k) = -nd
      k = k + 1
      nabor(k) = 0
      nb = nb + 1
    END IF
    
!   Increment NL and print X, Y, and NABOR.
    
    inc = (k-1)/8 + 2
    nl = nl + inc
    IF (nl > nlmax) THEN
      WRITE (lun,106)
      nl = inc
    END IF
    WRITE (lun,104) node, x(node), y(node), nabor(1:k)
    IF (k /= 8) WRITE (lun,105)
  END DO
END IF

! Print NB, NA, and NT (boundary nodes, arcs, and triangles).

nt = 2*nn - nb - 2
na = nt + nn - 1
IF (nl > nlmax-6) WRITE (lun,106)
WRITE (lun,107) nb, na, nt

! Print NCC and LCC.

5 WRITE (lun,108) ncc
IF (ncc > 0) WRITE (lun,109) lcc(1:ncc)
RETURN

! Print formats:

100 FORMAT (//t27, 'Adjacency Sets,    N = ', i5//)
101 FORMAT (' Node', t38, 'Neighbours of Node'//)
102 FORMAT (' Node     X(Node)        Y(Node)', t53, 'Neighbours of Node'//)
103 FORMAT (' ', i4, '     ', 14I5 / (t11, 14I5))
104 FORMAT (' ', i4, 2g15.6, '     ', 8I5 / (t41, 8I5))
105 FORMAT (' ')
106 FORMAT (//)
107 FORMAT (/' NB = ', i4, ' Boundary Nodes     NA = ', i5,  &
             ' Arcs     NT = ', i5, ' Triangles')
108 FORMAT (/' NCC =', i3, ' Constraint Curves')
109 FORMAT (t11, 14I5)
110 FORMAT (t12, '*** N is outside its valid', ' range ***')
END SUBROUTINE trprnt



FUNCTION frwrd(xa, ya, xb, yb, xc, yc) RESULT(fn_val)
REAL, INTENT(IN)  :: xa, ya, xb, yb, xc, yc
LOGICAL           :: fn_val

fn_val = (xb-xa)*(xc-xa) + (yb-ya)*(yc-ya) >= 0.
RETURN
END FUNCTION frwrd

END MODULE TRIPACK



PROGRAM Test_TriPack

!        TRITEST:  Portable Test Driver for TRIPACK
!                        07/02/98

!   This driver tests software package TRIPACK for constructing a constrained
! Delaunay triangulation of a set of points in the plane.  All modules other
! than TRMSHR are tested unless an error is encountered, in which case the
! program terminates immediately.

!   By default, tests are performed on a simple data set consisting of 12
! nodes whose convex hull covers the unit square.  The data set includes a
! single constraint region consisting of four nodes forming a smaller square
! at the center of the unit square.  However, by enabling the READ statements
! below (C# in the first two columns), testing may be performed on an
! arbitrary set of up to NMAX nodes with up to NCMAX constraint curves.
! (Refer to the PARAMETER statements below.)  A data set consists of the
! following sequence of records:

!    N = Number of nodes (format I4) -- 3 to NMAX.
!    NCC = Number of constraint curves (format I4) -- 0 to NCMAX.
!    (LCC(I), I = 1,NCC) = Indexes of the first node in each constraint curve
!             (format I4).  1 <= LCC(1) and, for I > 1, LCC(I-1) + 3 <= LCC(I)
!             <= N-2. (Each constraint curve has at least three nodes.)
!    (X(I),Y(I), I = 1,N) = Nodal coordinates with non-constraint nodes
!                followed by the NCC sequences of constraint nodes (format
!                2F13.8).

!   The I/O units may be changed by altering LIN (input) and
! LOUT (output) in the DATA statement below.

!   This driver must be linked to TRIPACK.

! Code converted using TO_F90 by Alan Miller
! Date: 2000-05-25  Time: 18:21:14

USE TriPack
IMPLICIT NONE

CHARACTER (LEN=80)  :: title
INTEGER             :: ier, io1, io2, k, ksum, lnew, lw, n0, na, nb, nn, nt
LOGICAL             :: numbr, prntx
REAL                :: a, armax, dsq, wx1, wx2, wy1, wy2

INTEGER, PARAMETER  :: nmax=100, ncmax=5, ntmx=2*nmax, n6=6*nmax,  &
                       lwk=2*nmax, nrow=9

! Array storage:

INTEGER  :: lcc(ncmax), lct(ncmax), lend(nmax), list(n6),  &
            lptr(n6), ltri(nrow,ntmx), nodes(lwk), iwk(2,nmax)
REAL     :: ds(nmax), x(nmax), y(nmax)

! Tolerance for TRMTST and NEARND:  upper bound on squared distances.

REAL, PARAMETER  :: tol = 1.0e-2

! Plot size for the triangulation plot.

REAL, PARAMETER  :: pltsiz = 7.5

! Logical unit numbers for I/O:

INTEGER, PARAMETER  :: lin = 1,  lout = 2,  lplt = 3

! Default data set:

INTEGER  :: n = 12, ncc = 1

x(1:12) = (/ 0., 1., .5, .15, .85, .5, 0., 1., .35, .65, .65, .35 /)
y(1:12) = (/ 0., 0., .15, .5, .5, .85, 1., 1., .35, .35, .65, .65 /)

OPEN (lout, FILE='RES')
OPEN (lplt, FILE='RES.eps')

! Store a plot title.  It must be enclosed in parentheses.

title = '(Triangulation created by TRITEST)'

! *** Read triangulation parameters -- N, NCC, LCC, X, Y.

!#    OPEN (LIN,FILE='tritest.dat',STATUS='OLD')
!#    READ (LIN,100,ERR=30) N, NCC
IF (n < 3 .OR. n > nmax .OR. ncc < 0 .OR. ncc > ncmax) GO TO 31
!#    IF (NCC .GT. 0) READ (LIN,100,ERR=30)
!#   .                     (LCC(K), K = 1,NCC)
!#    READ (LIN,110,ERR=30) (X(K),Y(K), K = 1,N)
!#100 FORMAT (I4)
!#110 FORMAT (2F13.8)

lcc(1) = 9

! Print a heading.

WRITE (lout,400) n

! *** Create the Delaunay triangulation (TRMESH), and test for errors
!     (refer to TRMTST below).  NODES and DS are used as work space.

CALL trmesh (n,x,y, list,lptr,lend,lnew,nodes, nodes(n+1:),ds,ier)
IF (ier == -2) THEN
  WRITE (lout,210)
  STOP
ELSE IF (ier == -4) THEN
  WRITE (lout,212)
  STOP
ELSE IF (ier > 0) THEN
  WRITE (lout,214)
  STOP
END IF
CALL trmtst (n,x,y,list,lptr,lend,lnew,tol, lout, armax,ier)
WRITE (lout,410) armax
IF (ier > 0) STOP

! *** Add the constraint curves (ADDCST).  Note that edges and triangles are
!     not removed from constraint regions.
!     ADDCST forces the inclusion of triangulation edges connecting the
!     sequences of constraint nodes.  If it is necessary to alter the
!     triangulation, the empty circumcircle property is no longer satisfied.

lw = lwk
CALL addcst (ncc,lcc,n,x,y, lw,list,lptr, lend, ier)
IF (ier /= 0) THEN
  WRITE (lout,220) ier
  STOP
END IF
IF (lw == 0) WRITE (lout,430)

! *** Test TRPRNT, TRLIST, and TRLPRT, and TRPLOT.

prntx = .true.
CALL trprnt (ncc,lcc,n,x,y,list,lptr,lend,lout,prntx)
CALL trlist (ncc,lcc,n,list,lptr,lend,nrow, nt,ltri, lct,ier)
CALL trlprt (ncc,lct,n,x,y,nrow,nt,ltri,lout,prntx)

!   Set the plot window [WX1,WX2] X [WY1,WY2] to the
!     smallest rectangle that contains the nodes.
!     NUMBR = TRUE iff nodal indexes are to be displayed.

wx1 = x(1)
wx2 = wx1
wy1 = y(1)
wy2 = wy1
DO  k = 2,n
  IF (x(k) < wx1) wx1 = x(k)
  IF (x(k) > wx2) wx2 = x(k)
  IF (y(k) < wy1) wy1 = y(k)
  IF (y(k) > wy2) wy2 = y(k)
END DO
numbr = n <= 200
CALL trplot (lplt,pltsiz,wx1,wx2,wy1,wy2,ncc,lcc,n,  &
             x,y,list,lptr,lend,title,numbr, ier)
IF (ier == 0) THEN
  WRITE (lout,470)
ELSE
  WRITE (lout,270) ier
END IF

! *** Test BNODES and AREAP.

CALL bnodes (n,list,lptr,lend, nodes,nb,na,nt)
a = areap(x,y,nb,nodes)
WRITE (lout,420) nb, na, nt, a

! *** Test GETNP by ordering the nodes on distance from N0
!                and verifying the ordering.

n0 = n/2
nodes(1) = n0
ds(1) = 0.
ksum = n0
DO  k = 2,n
  CALL getnp (ncc,lcc,n,x,y,list,lptr,lend, k, nodes,ds, ier)
  IF (ier /= 0 .OR. ds(k) < ds(k-1)) THEN
    WRITE (lout,230)
    STOP
  END IF
  ksum = ksum + nodes(k)
END DO

!   Test for all nodal indexes included in NODES.

IF (ksum /= (n*(n+1))/2) THEN
  WRITE (lout,230)
  STOP
END IF

! *** Test NEARND by verifying that the nearest node to K is
!                 node K for K = 1 to N.

DO  k = 1,n
  n0 = nearnd (x(k),y(k),1,n,x,y,list,lptr,lend, dsq)
  IF (n0 /= k  .OR.  dsq > tol) THEN
    WRITE (lout,240)
    STOP
  END IF
END DO

! *** Test DELARC by removing a boundary arc if possible.
!                 The first two nodes define a boundary arc
!                 in the default data set.

io1 = 1
io2 = 2
CALL delarc (n,io1,io2, list,lptr,lend,lnew, ier)
IF (ier == 1 .OR. ier == 4) THEN
  WRITE (lout,250) ier
  STOP
END IF
IF (ier /= 0) WRITE (lout,440)

!   Recreate the triangulation without constraints.

CALL trmesh (n,x,y, list,lptr,lend,lnew,nodes, nodes(n+1:),ds,ier)
ncc = 0

! *** Test DELNOD by removing nodes 4 to N (in reverse order).

IF (n <= 3) THEN
  WRITE (lout,450)
ELSE
  nn = n
  4 lw = lwk/2
  CALL delnod (nn,ncc, lcc,nn,x,y,list,lptr,lend, lnew,lw,iwk, ier)
  IF (ier /= 0) THEN
    WRITE (lout,260) ier
    STOP
  END IF
  IF (nn > 3) GO TO 4
END IF

! Successful test.

WRITE (lout,460)
STOP

! Error reading the data set.

!# 30 WRITE (*,200)
STOP

! Invalid value of N or NCC.

31 WRITE (*,205) n, ncc
STOP

! Error message formats:

!#200 FORMAT (//5X,'*** Input data set invalid ***'/)
205 FORMAT (//'     *** N or NCC is outside its valid ',  &
            'range:  N =', i5, ', NCC = ', i4, ' ***'/)
210 FORMAT (//'     *** Error in TRMESH:  the first three ',  &
            'nodes are collinear ***'/)
212 FORMAT (//'     *** Error in TRMESH:  invalid ', 'triangulation ***'/)
214 FORMAT (//'     *** Error in TRMESH:  duplicate nodes ', 'encountered ***'/)
220 FORMAT (//'     *** Error in ADDCST:  IER = ',i1, ' ***'/)
230 FORMAT (//'     *** Error in GETNP ***'/)
240 FORMAT (//'     *** Error in NEARND ***'/)
250 FORMAT (//'     *** Error in DELARC:  IER = ',i1, ' ***'/)
260 FORMAT (//'     *** Error in DELNOD:  IER = ',i1, ' ***'/)
270 FORMAT (//'     *** Error in TRPLOT:  IER = ',i1, ' ***'/)

! Informative message formats:

400 FORMAT (///t23, 'TRIPACK Test:  N =', i5///)
410 FORMAT ('     Maximum triangle aspect ratio = ', e10.3//)
420 FORMAT (/'     Output from BNODES and AREAP'//  &
            '     BNODES:  ', i4, ' boundary nodes,  ', i5, ' edges,  ', i5, &
            ' triangles'/  &
            '     AREAP:  area of convex hull = ', e15.8//)
430 FORMAT ('     Subroutine EDGE not tested:'/  &
            '       No edges were swapped by ADDCST'//)
440 FORMAT ('     Subroutine DELARC not tested:'/  &
            '       Nodes 1 and 2 do not form a ', 'removable boundary edge.'//)
450 FORMAT ('     Subroutine DELNOD not tested:'/  &
    '       N cannot be reduced below 3'//)
460 FORMAT ('     No triangulation errors encountered.'/)
470 FORMAT (/'     A triangulation plot file was ', 'successfully created.'/)

CONTAINS


SUBROUTINE trmtst (n,x,y,list,lptr,lend,lnew,tol, lun, armax,ier)

INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN)      :: x(:)
REAL, INTENT(IN)      :: y(:)
INTEGER, INTENT(IN)   :: list(:)
INTEGER, INTENT(IN)   :: lptr(:)
INTEGER, INTENT(IN)   :: lend(:)
INTEGER, INTENT(IN)   :: lnew
REAL, INTENT(IN)      :: tol
INTEGER, INTENT(IN)   :: lun
REAL, INTENT(OUT)     :: armax
INTEGER, INTENT(OUT)  :: ier

!***********************************************************

!                                               From TRPTEST
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                             (817) 565-2767
!                                                   09/01/91

!   This subroutine tests the validity of the data structure representing a
! Delaunay triangulation created by subroutine TRMESH.  The following
! properties are tested:

!   1)  Each interior node has at least three neighbors, and
!       each boundary node has at least two neighbors.

!   2)  abs(LIST(LP)) is a valid nodal index in the range 1 to N and
!       LIST(LP) > 0 unless LP = LEND(K) for some nodal index K.

!   3)  Each pointer LEND(K) for K = 1 to N and LPTR(LP) for
!       LP = 1 to LNEW-1 is a valid LIST index in the range 1 to LNEW-1.

!   4)  N .GE. NB .GE. 3, NT = 2*N-NB-2, and NA = 3*N-NB-3 = (LNEW-1)/2,
!       where NB, NT, and NA are the numbers of boundary nodes,
!       triangles, and arcs, respectively.

!   5)  Each circumcircle defined by the vertices of a triangle contains no
!       nodes in its interior.  This property distinguishes a Delaunay
!       triangulation from an arbitrary triangulation of the nodes.

! Note that no test is made for the property that a triangulation covers
! the convex hull of the nodes.

! On input:

!       N = Number of nodes.  N >= 3.

!       X,Y = Arrays of length N containing the nodal coordinates.

!       LIST,LPTR,LEND = Data structure containing the triangulation.
!                        Refer to subroutine TRMESH.

!       TOL = Nonnegative tolerance to allow for floating-point errors in the
!             circumcircle test.  An error situation is defined as
!             (R**2 - D**2)/R**2 > TOL, where R is the radius of a
!             circumcircle and D is the distance from the circumcenter to the
!             nearest node.  A reasonable value for TOL is 10*EPS, where EPS
!             is the machine precision.  The test is effectively bypassed by
!             making TOL large.  If TOL < 0, the tolerance is taken to be 0.

!       LUN = Logical unit number for printing error messages.
!             If LUN < 0 or LUN > 99, no messages are printed.

! Input parameters are not altered by this routine.

! On output:

!       ARMAX = Maximum aspect ratio (radius of inscribed circle divided by
!               circumradius) of a triangle in the triangulation unless
!               IER > 0.

!       IER = Error indicator:
!             IER = -1 if one or more null triangles (area = 0) are present
!                      but no (other) errors were encountered.  A null
!                      triangle is an error only if it occurs in the the
!                      interior.
!             IER = 0 if no errors or null triangles were encountered.
!             IER = 1 if a node has too few neighbors.
!             IER = 2 if a LIST entry is outside its valid range.
!             IER = 3 if a LPTR or LEND entry is outside its valid range.
!             IER = 4 if the triangulation parameters (N, NB, NT, NA, and
!                     LNEW) are inconsistent (or N < 3 or LNEW is invalid).
!             IER = 5 if a triangle contains a node interior to its
!                     circumcircle.

! Module required by TRMTST:  CIRCUM

! Intrinsic function called by TRMTST:  MAX, ABS

!***********************************************************

INTEGER :: k, lp, lpl, lpmax, lpn, n1, n2, n3, na, nb,  &
           nfail, nn, nnb, nt, null
LOGICAL :: ratio, rite
REAL    :: ar, cr, cx, cy, rs, rtol, sa

! Store local variables, test for errors in input, and initialize counts.

nn = n
lpmax = lnew - 1
rtol = tol
IF (rtol < 0.) rtol = 0.
ratio = .true.
armax = 0.
rite = lun >= 0 .AND. lun <= 99
IF (nn < 3) GO TO 14
nb = 0
nt = 0
null = 0
nfail = 0

! Loop on triangles (N1,N2,N3) such that N2 and N3 index adjacent neighbors of
!   N1 and are both larger than N1 (each triangle is associated with its
!   smallest index).
!   NNB is the neighbor count for N1.

DO  n1 = 1,nn
  nnb = 0
  lpl = lend(n1)
  IF (lpl < 1 .OR. lpl > lpmax) THEN
    lp = lpl
    GO TO 13
  END IF
  lp = lpl
  
!   Loop on neighbors of N1.
  
  1 lp = lptr(lp)
  nnb = nnb + 1
  IF (lp < 1 .OR. lp > lpmax) GO TO 13
  n2 = list(lp)
  IF (n2 < 0) THEN
    IF (lp /= lpl) GO TO 12
    IF (n2 == 0 .OR. -n2 > nn) GO TO 12
    nb = nb + 1
    GO TO 4
  END IF
  IF (n2 < 1 .OR. n2 > nn) GO TO 12
  lpn = lptr(lp)
  n3 = ABS(list(lpn))
  IF (n2 < n1 .OR. n3 < n1) GO TO 4
  nt = nt + 1
  
!   Compute the coordinates of the circumcenter of (N1,N2,N3).
  
  CALL circum (x(n1),y(n1),x(n2),y(n2),x(n3),y(n3), ratio, cx,cy,cr,sa,ar)
  IF (sa == 0.) THEN
    null = null + 1
    GO TO 4
  END IF
  armax = MAX(armax,ar)
  
!   Test for nodes within the circumcircle.
  
  rs = cr*cr*(1.-rtol)
  DO  k = 1,nn
    IF (k == n1 .OR. k == n2 .OR. k == n3) CYCLE
    IF ((cx-x(k))**2 + (cy-y(k))**2 < rs) GO TO 3
  END DO
  GO TO 4
  
!   Node K is interior to the circumcircle of (N1,N2,N3).
  
  3 nfail = nfail + 1
  
!   Bottom of loop on neighbors.
  
  4 IF (lp /= lpl) GO TO 1
  IF (nnb < 2 .OR. (nnb == 2 .AND. list(lpl) > 0)) GO TO 11
END DO

! Test parameters for consistency and check for NFAIL = 0.

na = lpmax/2
IF (nb < 3 .OR. nt /= 2*nn-nb-2 .OR. na /= 3*nn-nb-3) GO TO 14
IF (nfail /= 0) GO TO 15

! No errors were encountered.

ier = 0
IF (null == 0) RETURN
ier = -1
IF (rite) WRITE (lun,100) null
100 FORMAT (//'     *** TRMTST -- ',i5,' NULL TRIANGLES ARE PRESENT'/  &
            t20, '(NULL TRIANGLES ON THE BOUNDARY ARE UNAVOIDABLE) ***'//)
RETURN

! Node N1 has fewer than three neighbors.

11 ier = 1
IF (rite) WRITE (lun,110) n1, nnb
110 FORMAT (//'     *** TRMTST -- NODE ',i5, ' HAS ONLY ',i5,' NEIGHBORS ***'/)
RETURN

! N2 = LIST(LP) is outside its valid range.

12 ier = 2
IF (rite) WRITE (lun,120) n2, lp, n1
120 FORMAT (//'     *** TRMTST -- LIST(LP) =', i5, ', FOR LP =', i5, ','/  &
            t20, 'IS NOT A VALID NEIGHBOR OF ', i5, ' ***'/)
RETURN

! LIST pointer LP is outside its valid range.

13 ier = 3
IF (rite) WRITE (lun,130) lp, lnew, n1
130 FORMAT (//'     *** TRMTST -- LP =',i5,' IS NOT IN THE',  &
    ' RANGE 1 TO LNEW-1 FOR LNEW = ',i5/ 19X,'LP POINTS TO A NEIGHBOR OF ',i5,  &
    ' ***'/)
RETURN

! Inconsistent triangulation parameters encountered.

14 ier = 4
IF (rite) WRITE (lun,140) n, lnew, nb, nt, na
140 FORMAT (//'     *** TRMTST -- INCONSISTENT PARAMETERS',  &
    ' ***'/19X,'N = ',i5,' NODES',12X,'LNEW =',i5/  &
    19X,'NB = ',i5,' BOUNDARY NODES'/ 19X,'NT = ',i5,' TRIANGLES'/  &
    19X,'NA = ',i5,' ARCS'/)
RETURN

! Circumcircle test failure.

15 ier = 5
IF (rite) WRITE (lun,150) nfail
150 FORMAT (//'     *** TRMTST -- ',i5,' CIRCUMCIRCLES ',  &
    'CONTAIN NODES IN THEIR INTERIORS ***'/)
RETURN
END SUBROUTINE trmtst

END PROGRAM Test_TriPack


!--------------------------------------------------------------------------

! The output file `RES' follows



                      TRIPACK Test:  N =   12



     Maximum triangle aspect ratio =  0.480E+00




                          Adjacency Sets,    N =    12


 Node     X(Node)        Y(Node)                    Neighbours of Node


    1    0.00000        0.00000             2    3    9    4    7    0
 
    2    1.00000        0.00000             8    5   10    3    1    0
 
    3   0.500000       0.150000            10    9    1    2
 
    4   0.150000       0.500000             7    1    9   12
 
    5   0.850000       0.500000             8   11   10    2
 
    6   0.500000       0.850000             8    7   12   11
 
    7    0.00000        1.00000             1    4   12    6    8    0
 
    8    1.00000        1.00000             7    6   11    5    2    0
 
    9   0.350000       0.350000            10   11   12    4    1    3
 
   10   0.650000       0.350000             2    5   11    9    3
 
   11   0.650000       0.650000            12    9   10    5    8    6
 
   12   0.350000       0.650000             9   11    6    7    4
 

 NB =    4 Boundary Nodes     NA =    29 Arcs     NT =    18 Triangles

 NCC =  1 Constraint Curves
              9



                        TRIPACK (TRLIST) Output


                Node       X(Node)          Y(Node)


                   1     0.000000E+00     0.000000E+00
                   2     0.100000E+01     0.000000E+00
                   3     0.500000E+00     0.150000E+00
                   4     0.150000E+00     0.500000E+00
                   5     0.850000E+00     0.500000E+00
                   6     0.500000E+00     0.850000E+00
                   7     0.000000E+00     0.100000E+01
                   8     0.100000E+01     0.100000E+01
                   9     0.350000E+00     0.350000E+00
                  10     0.650000E+00     0.350000E+00
                  11     0.650000E+00     0.650000E+00
                  12     0.350000E+00     0.650000E+00


 Triangle        Vertices            Neighbors              Arcs
    KT       N1     N2     N3    KT1    KT2    KT3    KA1    KA2    KA3

     1        1      2      3      7      2      0      8      2      1
     2        1      3      9      8      3      1     10      3      2
     3        1      9      4      9      4      2     12      5      3
     4        1      4      7     10      0      3     13      4      5
     5        2      8      5     11      6      0     15      7      6
     6        2      5     10     12      7      5     16      9      7
     7        2     10      3      8      1      6     11      8      9
     8        3     10      9     17      2      7     26     10     11
     9        4      9     12     18     10      3     28     14     12
    10        4     12      7     14      4      9     19     13     14
    11        5      8     11     16     12      5     22     17     15
    12        5     11     10     17      6     11     25     16     17
    13        6      8      7      0     14     16     18     20     23
    14        6      7     12     10     15     13     19     21     20
    15        6     12     11     18     16     14     27     24     21
    16        6     11      8     11     13     15     22     23     24
    17        9     10     11     12     18      8     25     29     26
    18        9     11     12     15      9     17     27     28     29

 NB =    4 Boundary Nodes     NA =    29 Arcs     NT =    18 Triangles

 NCC =  1 Constraint Curves
             17

     A triangulation plot file was successfully created.


     Output from BNODES and AREAP

     BNODES:     4 boundary nodes,     29 edges,     18 triangles
     AREAP:  area of convex hull =  0.10000000E+01


     No triangulation errors encountered.


!--------------------------------------------------------------------------

! The output file `RES.EPS' follows


%!PS-Adobe-3.0 EPSF-3.0
%%BoundingBox:  36 126 576 666
%%Title:  Triangulation
%%Creator:  TRIPACK
%%EndComments
    2.000000 setlinewidth
  69 159 moveto
  69 633 lineto
 543 633 lineto
 543 159 lineto
closepath
stroke
   69.000000  159.000000 translate
  475.000000  475.000000 scale
    0.002105 setlinewidth
gsave
    0.000000    0.000000 moveto
    1.000000    0.000000 lineto
    1.000000    1.000000 lineto
    0.000000    1.000000 lineto
closepath clip newpath
    0.000000    0.000000 moveto    1.000000    0.000000 lineto
    0.000000    0.000000 moveto    0.500000    0.150000 lineto
    0.000000    0.000000 moveto    0.350000    0.350000 lineto
    0.000000    0.000000 moveto    0.150000    0.500000 lineto
    0.000000    0.000000 moveto    0.000000    1.000000 lineto
    1.000000    0.000000 moveto    1.000000    1.000000 lineto
    1.000000    0.000000 moveto    0.850000    0.500000 lineto
    1.000000    0.000000 moveto    0.650000    0.350000 lineto
    1.000000    0.000000 moveto    0.500000    0.150000 lineto
    0.500000    0.150000 moveto    0.650000    0.350000 lineto
    0.500000    0.150000 moveto    0.350000    0.350000 lineto
    0.150000    0.500000 moveto    0.000000    1.000000 lineto
    0.150000    0.500000 moveto    0.350000    0.350000 lineto
    0.150000    0.500000 moveto    0.350000    0.650000 lineto
    0.850000    0.500000 moveto    1.000000    1.000000 lineto
    0.850000    0.500000 moveto    0.650000    0.650000 lineto
    0.850000    0.500000 moveto    0.650000    0.350000 lineto
    0.500000    0.850000 moveto    1.000000    1.000000 lineto
    0.500000    0.850000 moveto    0.000000    1.000000 lineto
    0.500000    0.850000 moveto    0.350000    0.650000 lineto
    0.500000    0.850000 moveto    0.650000    0.650000 lineto
    0.000000    1.000000 moveto    0.350000    0.650000 lineto
    0.000000    1.000000 moveto    1.000000    1.000000 lineto
    1.000000    1.000000 moveto    0.650000    0.650000 lineto
stroke
[    0.008421] 0 setdash
    0.350000    0.350000 moveto    0.650000    0.350000 lineto
    0.350000    0.350000 moveto    0.650000    0.650000 lineto
    0.350000    0.350000 moveto    0.350000    0.650000 lineto
    0.650000    0.350000 moveto    0.650000    0.650000 lineto
    0.650000    0.650000 moveto    0.350000    0.650000 lineto
stroke
grestore
/Helvetica findfont
    0.021053 scalefont setfont
    0.000000    0.000000 moveto
(  1) show
    1.000000    0.000000 moveto
(  2) show
    0.500000    0.150000 moveto
(  3) show
    0.150000    0.500000 moveto
(  4) show
    0.850000    0.500000 moveto
(  5) show
    0.500000    0.850000 moveto
(  6) show
    0.000000    1.000000 moveto
(  7) show
    1.000000    1.000000 moveto
(  8) show
    0.350000    0.350000 moveto
(  9) show
    0.650000    0.350000 moveto
( 10) show
    0.650000    0.650000 moveto
( 11) show
    0.350000    0.650000 moveto
( 12) show
/Helvetica findfont
    0.033684 scalefont setfont
(Triangulation created by TRITEST)                                              
  stringwidth pop 2 div neg     0.500000 add     1.101053 moveto
(Triangulation created by TRITEST)                                              
  show
    0.000000   -0.105263 moveto
(Window:   WX1 = 0.000E+00,   WX2 = 0.100E+01) show
(Window:  ) stringwidth pop     0.000000 add   -0.172632 moveto
( WY1 = 0.000E+00,   WY2 = 0.100E+01) show
stroke
showpage
%%EOF
