MODULE coloring
! A module from MINPACK-2, used by DTRON

IMPLICIT NONE

INTEGER, PARAMETER, PRIVATE :: dp = SELECTED_REAL_KIND(14, 60)


CONTAINS


SUBROUTINE dssm(n, npairs, indrow, indcol, method, listp, ngrp,  &
                maxgrp, mingrp, info, ipntr, jpntr)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-30  Time: 12:40:40
! Latest revision - 21 November 1999

INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN)      :: npairs
INTEGER, INTENT(IN OUT)  :: indrow(:)    ! indrow(npairs)
INTEGER, INTENT(IN OUT)  :: indcol(:)    ! indcol(npairs)
INTEGER, INTENT(IN)      :: method
INTEGER, INTENT(OUT)     :: listp(:)
INTEGER, INTENT(OUT)     :: ngrp(:)
INTEGER, INTENT(OUT)     :: maxgrp
INTEGER, INTENT(OUT)     :: mingrp
INTEGER, INTENT(OUT)     :: info
INTEGER, INTENT(IN OUT)  :: ipntr(:)    ! ipntr(n+1)
INTEGER, INTENT(IN OUT)  :: jpntr(:)    ! jpntr(n+1)

!  **********

!  subroutine dssm

!  Given the sparsity pattern of a symmetric matrix A of order n,
!  this subroutine determines a symmetric permutation of A and a
!  partition of the columns of A consistent with the determination
!  of A by a lower triangular substitution method.

!  The sparsity pattern of the matrix A is specified by the arrays indrow and
!  indcol.  On input the indices for the non-zero elements in the lower
!  triangular part of A are

!        (indrow(k),indcol(k)), k = 1,2,...,npairs.

!  The (indrow(k),indcol(k)) pairs may be specified in any order.
!  Duplicate input pairs are permitted, but the subroutine eliminates them.
!  The subroutine requires that all the diagonal elements be part of the
!  sparsity pattern and replaces any pair (indrow(k),indcol(k)) where indrow(k)
!  is less than indcol(k) by the pair (indcol(k),indrow(k)).

!  The direct method (method = 1) first determines a partition of the columns
!  of A such that two columns in a group have a non-zero element in row k only
!  if column k is in an earlier group.  Using this partition, the subroutine
!  then computes a symmetric permutation of A consistent with the determination
!  of A by a lower triangular substitution method.

!  The indirect method first computes a symmetric permutation of A which
!  minimizes the maximum number of non-zero elements in any row of L, where L
!  is the lower triangular part of the permuted matrix.  The subroutine then
!  partitions the columns of L into groups such that columns of L in a group
!  do not have a non-zero in the same row position.

!  The subroutine statement is

!    subroutine dssm(n, npairs, indrow, indcol, method, listp, ngrp,
!                    maxgrp, mingrp, info, ipntr, jpntr)

!  where

!    n is a positive integer input variable set to the order of A.

!    npairs is a positive integer input variable set to the number of
!      (indrow,indcol) pairs used to describe the sparsity pattern of A.

!    indrow is an integer array of length npairs.  On input indrow must contain
!      the row indices of the non-zero elements in the lower triangular part
!      of A.  On output indrow is permuted so that the corresponding column
!      indices are in non-decreasing order.  The column indices can be
!      recovered from the array jpntr.

!    indcol is an integer array of length npairs.  On input indcol must contain
!      the column indices of the non-zero elements in the lower triangular part
!      of A.  On output indcol is permuted so that the corresponding row
!      indices are in non-decreasing order.  The row indices can be recovered
!      from the array ipntr.

!    method is an integer input variable.  If method = 1, the direct method
!      is used to determine the partition and symmetric permutation.
!      Otherwise, the indirect method is used to determine the symmetric
!      permutation and partition.

!    listp is an integer output array of length n which specifies
!      the symmetric permutation of the matrix A.  Element (i,j)
!      of A is the (listp(i),listp(j)) element of the permuted matrix.

!    ngrp is an integer output array of length n which specifies the
!      partition of the columns of A. Column j belongs to group ngrp(j).

!    maxgrp is an integer output variable which specifies the
!      number of groups in the partition of the columns of A.

!    mingrp is an integer output variable which specifies a lower bound for
!      the number of groups in any partition of the columns of A consistent
!      with the determination of A by a lower triangular substitution method.

!    info is an integer output variable set as follows. For
!      normal termination info = 1. If n or npairs is not
!      positive or liwa is less than 6*n, then info = 0. If the
!      k-th element of indrow or the k-th element of indcol is
!      not an integer between 1 and n, or if the k-th diagonal
!      element is not in the sparsity pattern, then info = -k.

!    ipntr is an integer output array of length n + 1 which specifies the
!      locations of the column indices in indcol.
!      The column indices for row i are

!            indcol(k), k = ipntr(i),...,ipntr(i+1)-1.

!      Note that ipntr(n+1)-1 is then the number of non-zero
!      elements in the lower triangular part of the matrix A.

!    jpntr is an integer output array of length n + 1 which specifies the
!      locations of the row indices in indrow.
!      The row indices for column j are

!            indrow(k), k = jpntr(j),...,jpntr(j+1)-1.

!      Note that jpntr(n+1)-1 is then the number of non-zero
!      elements in the lower triangular part of the matrix A.

!  Subprograms called

!    MINPACK-supplied ... degr,ido,idog,numsrt,sdpt,seq,setr,slo,slog,srtdat

!    FORTRAN-supplied ... max,min

!  Argonne National Laboratory. MINPACK Project. December 1984.
!  Thomas F. Coleman, Burton S. Garbow, Jorge J. More'

!  **********

INTEGER :: i, ir, j, jp, k, maxid, maxvd, maxclq, nnz, numgrp
INTEGER :: iwa(3*n)

!     Check the input data.

info = 0
IF (n < 1 .OR. npairs < 1) RETURN
iwa(1:n) = 0
DO  k = 1, npairs
  info = -k
  IF (indrow(k) < 1 .OR. indrow(k) > n .OR.  &
      indcol(k) < 1 .OR. indcol(k) > n) RETURN
  IF (indrow(k) == indcol(k)) iwa(indrow(k)) = 1
END DO
DO  k = 1, n
  info = -k
  IF (iwa(k) /= 1) RETURN
END DO
info = 1

!     Generate the sparsity pattern for the lower triangular part of A.

DO  k = 1, npairs
  i = indrow(k)
  j = indcol(k)
  indrow(k) = MAX(i,j)
  indcol(k) = MIN(i,j)
END DO

!     Sort the data structure by columns.

CALL srtdat(n, npairs, indrow, indcol, jpntr)

!     Compress the data and determine the number of non-zero elements
!     in the lower triangular part of A.

iwa(1:n) = 0
nnz = 0
DO  j = 1, n
  k = nnz
  DO  jp = jpntr(j), jpntr(j+1)-1
    ir = indrow(jp)
    IF (iwa(ir) /= j) THEN
      nnz = nnz + 1
      indrow(nnz) = ir
      iwa(ir) = j
    END IF
  END DO
  jpntr(j) = k + 1
END DO
jpntr(n+1) = nnz + 1

!     Extend the data structure to rows.

CALL setr(n, n, indrow, jpntr, indcol, ipntr)

!     Determine the smallest-last ordering of the vertices of the adjacency
!     graph of A, and from it determine a lower bound for the number of groups.

CALL slog(n, indrow, jpntr, indcol, ipntr, iwa, maxclq, maxvd)
mingrp =  1 + maxvd

!     Use the selected method.

IF (method == 1) THEN
  
!        Direct method. Determine a partition of the columns
!        of A by the Powell-Toint method.
  
  CALL sdpt(n, indrow, jpntr, indcol, ipntr, ngrp, maxgrp)
  
!        Define a symmetric permutation of A according to the
!        ordering of the column group numbers in the partition.
  
  CALL numsrt(n, maxgrp, ngrp, 1, iwa, iwa(2*n+1:), iwa(n+1:))
  DO  i = 1, n
    listp(iwa(i)) = i
  END DO
ELSE
  
!        Indirect method. Determine the incidence degree ordering of the
!        vertices of the adjacency graph of A and, together with the
!        smallest-last ordering, define a symmetric permutation of A.
  
  CALL idog(n, indrow, jpntr, indcol, ipntr, listp, maxclq, maxid)
  IF (maxid > maxvd) listp(1:n) = iwa(1:n)
  
!        Generate the sparsity pattern for the lower
!        triangular part L of the permuted matrix.
  
  DO  j = 1, n
    DO  jp = jpntr(j), jpntr(j+1)-1
      i = indrow(jp)
      indrow(jp) = MAX(listp(i),listp(j))
      indcol(jp) = MIN(listp(i),listp(j))
    END DO
  END DO
  
!        Sort the data structure by columns.
  
  CALL srtdat(n, nnz, indrow, indcol, jpntr)
  
!        Extend the data structure to rows.
  
  CALL setr(n, n, indrow, jpntr, indcol, ipntr)
  
!        Determine the degree sequence for the intersection
!        graph of the columns of L.
  
  CALL degr(n, indrow, jpntr, indcol, ipntr, iwa(2*n+1:))
  
!        Color the intersection graph of the columns of L
!        with the smallest-last (SL) ordering.
  
  CALL slo(n, indrow, jpntr, indcol, ipntr, iwa(2*n+1:), iwa(n+1:), maxclq)
  CALL seq(n, indrow, jpntr, indcol, ipntr, iwa(n+1:), iwa, maxgrp)
  DO  j = 1, n
    ngrp(j) = iwa(listp(j))
  END DO
  
!        Exit if the smallest-last ordering is optimal.
  
  IF (maxgrp /= maxclq) THEN
  
!        Color the intersection graph of the columns of L
!        with the incidence degree (ID) ordering.

    CALL ido(n, n, indrow, jpntr, indcol, ipntr, iwa(2*n+1:), iwa(n+1:), maxclq)
    CALL seq(n, indrow, jpntr, indcol, ipntr, iwa(n+1:), iwa, numgrp)

!        Retain the better of the two orderings.
  
    IF (numgrp < maxgrp) THEN
      maxgrp = numgrp
      DO  j = 1, n
        ngrp(j) = iwa(listp(j))
      END DO
    END IF
  END IF
  
!        Generate the sparsity pattern for the lower
!        triangular part of the original matrix.
  
  DO  j = 1, n
    iwa(listp(j)) = j
  END DO
  DO  j = 1, n
    DO  jp = jpntr(j), jpntr(j+1)-1
      i = indrow(jp)
      indrow(jp) = MAX(iwa(i),iwa(j))
      indcol(jp) = MIN(iwa(i),iwa(j))
    END DO
  END DO
  
!        Sort the data structure by columns.
  
  CALL srtdat(n, nnz, indrow, indcol, jpntr)
  
!        Extend the data structure to rows.
  
  CALL setr(n, n, indrow, jpntr, indcol, ipntr)
END IF
RETURN

!     Last card of subroutine dssm.

END SUBROUTINE dssm



SUBROUTINE degr(n, indrow, jpntr, indcol, ipntr, ndeg)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-30  Time: 12:48:53

INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: indrow(:)
INTEGER, INTENT(IN)   :: jpntr(:)    ! jpntr(n+1)
INTEGER, INTENT(IN)   :: indcol(:)
INTEGER, INTENT(IN)   :: ipntr(:)
INTEGER, INTENT(OUT)  :: ndeg(:)

!  **********

!  subroutine degr

!  Given the sparsity pattern of an m by n matrix A, this subroutine determines
!  the degree sequence for the intersection graph of the columns of A.

!  In graph-theory terminology, the intersection graph of the columns of A is
!  the loopless graph G with vertices a(j), j = 1,2,...,n where a(j) is the
!  j-th column of A and with edge (a(i),a(j)) if and only if columns i and j
!  have a non-zero in the same row position.

!  Note that the value of m is not needed by degr and is therefore not present
!  in the subroutine statement.

!  The subroutine statement is

!    subroutine degr(n, indrow, jpntr, indcol, ipntr, ndeg, iwa)

!  where

!    n is a positive integer input variable set to the number of columns of A.

!    indrow is an integer input array which contains the row
!      indices for the non-zeroes in the matrix A.

!    jpntr is an integer input array of length n + 1 which
!      specifies the locations of the row indices in indrow.
!      The row indices for column j are

!            indrow(k), k = jpntr(j),...,jpntr(j+1)-1.

!      Note that jpntr(n+1)-1 is then the number of non-zero
!      elements of the matrix A.

!    indcol is an integer input array which contains the
!      column indices for the non-zeroes in the matrix A.

!    ipntr is an integer input array of length m + 1 which
!      specifies the locations of the column indices in indcol.
!      The column indices for row i are

!            indcol(k), k = ipntr(i),...,ipntr(i+1)-1.

!      Note that ipntr(m+1)-1 is then the number of non-zero
!      elements of the matrix A.

!    ndeg is an integer output array of length n which specifies the degree
!      sequence.  The degree of the j-th column of A is ndeg(j).

!    iwa is an integer work array of length n.

!  Argonne National Laboratory. MINPACK Project. July 1983.
!  Thomas F. Coleman, Burton S. Garbow, Jorge J. More'

!  **********

INTEGER :: ic, ip, ir, jcol, jp, iwa(n)

!     Initialization block.

DO  jp = 1, n
  ndeg(jp) = 0
  iwa(jp) = 0
END DO

!     Compute the degree sequence by determining the contributions
!     to the degrees from the current(jcol) column and further
!     columns which have not yet been considered.

DO  jcol = 2, n
  iwa(jcol) = n
  
!        Determine all positions (ir,jcol) which correspond
!        to non-zeroes in the matrix.
  
  DO  jp = jpntr(jcol), jpntr(jcol+1)-1
    ir = indrow(jp)
    
!           For each row ir, determine all positions (ir,ic)
!           which correspond to non-zeroes in the matrix.
    
    DO  ip = ipntr(ir), ipntr(ir+1)-1
      ic = indcol(ip)
      
!              Array iwa marks columns which have contributed to
!              the degree count of column jcol. Update the degree
!              counts of these columns as well as column jcol.
      
      IF (iwa(ic) < jcol) THEN
        iwa(ic) = jcol
        ndeg(ic) = ndeg(ic) + 1
        ndeg(jcol) = ndeg(jcol) + 1
      END IF
    END DO
  END DO
END DO
RETURN

!     Last card of subroutine degr.

END SUBROUTINE degr



SUBROUTINE ido(m, n, indrow, jpntr, indcol, ipntr, ndeg, list, maxclq)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-30  Time: 12:48:58

INTEGER, INTENT(IN)      :: m
INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN)      :: indrow(:)
INTEGER, INTENT(IN)      :: jpntr(:)    ! jpntr(n+1)
INTEGER, INTENT(IN)      :: indcol(:)
INTEGER, INTENT(IN)      :: ipntr(:)    ! ipntr(m+1)
INTEGER, INTENT(IN)      :: ndeg(:)
INTEGER, INTENT(OUT)     :: list(:)
INTEGER, INTENT(OUT)     :: maxclq

!  **********

!  subroutine ido

!  Given the sparsity pattern of an m by n matrix A, this subroutine
!  determines an incidence-degree ordering of the columns of A.

!  The incidence-degree ordering is defined for the loopless
!  graph G with vertices a(j), j = 1,2,...,n where a(j) is the
!  j-th column of A and with edge (a(i),a(j)) if and only if
!  columns i and j have a non-zero in the same row position.

!  The incidence-degree ordering is determined recursively by letting
!  list(k), k = 1,...,n be a column with maximal incidence to the subgraph
!  spanned by the ordered columns.  Among all the columns of maximal incidence,
!  ido chooses a column of maximal degree.

!  The subroutine statement is

!    subroutine ido(m, n, indrow, jpntr, indcol, ipntr, ndeg, list,
!                   maxclq, iwa1, iwa2, iwa3, iwa4)

!  where

!    m is a positive integer input variable set to the number of rows of A.

!    n is a positive integer input variable set to the number of columns of A.

!    indrow is an integer input array which contains the row
!      indices for the non-zeroes in the matrix A.

!    jpntr is an integer input array of length n + 1 which
!      specifies the locations of the row indices in indrow.
!      The row indices for column j are

!            indrow(k), k = jpntr(j),...,jpntr(j+1)-1.

!      Note that jpntr(n+1)-1 is then the number of non-zero
!      elements of the matrix A.

!    indcol is an integer input array which contains the
!      column indices for the non-zeroes in the matrix A.

!    ipntr is an integer input array of length m + 1 which
!      specifies the locations of the column indices in indcol.
!      The column indices for row i are

!            indcol(k), k = ipntr(i),...,ipntr(i+1)-1.

!      Note that ipntr(m+1)-1 is then the number of non-zero
!      elements of the matrix A.

!    ndeg is an integer input array of length n which specifies the degree
!      sequence. The degree of the j-th column of A is ndeg(j).

!    list is an integer output array of length n which specifies
!      the incidence-degree ordering of the columns of A.  The j-th
!      column in this order is list(j).

!    maxclq is an integer output variable set to the size
!      of the largest clique found during the ordering.

!    iwa1, iwa2, iwa3, and iwa4 are integer work arrays of length n.

!  Subprograms called

!    MINPACK-supplied ... numsrt

!    FORTRAN-supplied ... max

!  Argonne National Laboratory. MINPACK Project. August 1984.
!  Thomas F. Coleman, Burton S. Garbow, Jorge J. More'

!  **********

INTEGER :: ic, ip, ir, jcol, jp, maxinc, maxlst, ncomp, numinc, numlst,  &
           numord, numwgt, iwa1(0:n-1), iwa2(n), iwa3(n), iwa4(n)

!     Sort the degree sequence.

CALL numsrt(n, n-1, ndeg, -1, iwa4, iwa2, iwa3)

!     Initialization block.

!     Create a doubly-linked list to access the incidences of the
!     columns.  The pointers for the linked list are as follows.

!     Each un-ordered column ic is in a list (the incidence list)
!     of columns with the same incidence.

!     iwa1(numinc) is the first column in the numinc list
!     unless iwa1(numinc) = 0.  In this case there are
!     no columns in the numinc list.

!     iwa2(ic) is the column before ic in the incidence list unless
!     iwa2(ic) = 0.
!     In this case ic is the first column in this incidence list.

!     iwa3(ic) is the column after ic in the incidence list unless iwa3(ic) = 0.
!     In this case ic is the last column in this incidence list.

!     If ic is an un-ordered column, then list(ic) is the incidence of ic to
!     the graph induced by the ordered columns.  If jcol is an ordered column,
!     then list(jcol) is the incidence-degree order of column jcol.

maxinc = 0
DO  jp = n, 1, -1
  ic = iwa4(jp)
  iwa1(n-jp) = 0
  iwa2(ic) = 0
  iwa3(ic) = iwa1(0)
  IF (iwa1(0) > 0) iwa2(iwa1(0)) = ic
  iwa1(0) = ic
  iwa4(jp) = 0
  list(jp) = 0
END DO

!     Determine the maximal search length for the list
!     of columns of maximal incidence.

maxlst = 0
DO  ir = 1, m
  maxlst = maxlst + (ipntr(ir+1) - ipntr(ir))**2
END DO
maxlst = maxlst/n
maxclq = 0
numord = 1

!     Beginning of iteration loop.

!        Choose a column jcol of maximal degree among the
!        columns of maximal incidence maxinc.

30 DO
  jp = iwa1(maxinc)
  IF (jp > 0) EXIT
  maxinc = maxinc - 1
END DO
numwgt = -1
DO  numlst = 1, maxlst
  IF (ndeg(jp) > numwgt) THEN
    numwgt = ndeg(jp)
    jcol = jp
  END IF
  jp = iwa3(jp)
  IF (jp <= 0) EXIT
END DO
list(jcol) = numord

!        Update the size of the largest clique found during the ordering.

IF (maxinc == 0) ncomp = 0
ncomp = ncomp + 1
IF (maxinc+1 == ncomp) maxclq = MAX(maxclq,ncomp)

!        Termination test.

numord = numord + 1
IF (numord > n) GO TO 100

!        Delete column jcol from the maxinc list.

IF (iwa2(jcol) == 0) THEN
  iwa1(maxinc) = iwa3(jcol)
ELSE
  iwa3(iwa2(jcol)) = iwa3(jcol)
END IF
IF (iwa3(jcol) > 0) iwa2(iwa3(jcol)) = iwa2(jcol)

!        Find all columns adjacent to column jcol.

iwa4(jcol) = n

!        Determine all positions (ir,jcol) which correspond
!        to non-zeroes in the matrix.

DO  jp = jpntr(jcol), jpntr(jcol+1)-1
  ir = indrow(jp)
  
!           For each row ir, determine all positions (ir,ic)
!           which correspond to non-zeroes in the matrix.
  
  DO  ip = ipntr(ir), ipntr(ir+1)-1
    ic = indcol(ip)
    
!              Array iwa4 marks columns which are adjacent to column jcol.
    
    IF (iwa4(ic) < numord) THEN
      iwa4(ic) = numord
      
!                 Update the pointers to the current incidence lists.
      
      numinc = list(ic)
      list(ic) = list(ic) + 1
      maxinc = MAX(maxinc,list(ic))
      
!                 Delete column ic from the numinc list.
      
      IF (iwa2(ic) == 0) THEN
        iwa1(numinc) = iwa3(ic)
      ELSE
        iwa3(iwa2(ic)) = iwa3(ic)
      END IF
      IF (iwa3(ic) > 0) iwa2(iwa3(ic)) = iwa2(ic)
      
!                 Add column ic to the numinc+1 list.
      
      iwa2(ic) = 0
      iwa3(ic) = iwa1(numinc+1)
      IF (iwa1(numinc+1) > 0) iwa2(iwa1(numinc+1)) = ic
      iwa1(numinc+1) = ic
    END IF
  END DO
END DO

!        End of iteration loop.

GO TO 30

!     Invert the array list.

100 DO  jcol = 1, n
  iwa2(list(jcol)) = jcol
END DO
list(1:n) = iwa2(1:n)
RETURN

!     Last card of subroutine ido.

END SUBROUTINE ido



SUBROUTINE idog(n, nghbrp, npntrp, nghbrs, npntrs, listp, maxclq, maxid)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-30  Time: 12:49:04

INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: nghbrp(:)
INTEGER, INTENT(IN)   :: npntrp(:)    ! npntrp(n+1)
INTEGER, INTENT(IN)   :: nghbrs(:)
INTEGER, INTENT(IN)   :: npntrs(:)    ! npntrs(n+1)
INTEGER, INTENT(OUT)  :: listp(:)
INTEGER, INTENT(OUT)  :: maxclq
INTEGER, INTENT(OUT)  :: maxid

!  **********

!  subroutine idog

!  Given a loopless graph G = (V,E), this subroutine determines
!  the incidence degree ordering of the vertices of G.

!  The incidence degree ordering is determined recursively by
!  letting list(k), k = 1,...,n be a vertex with maximal
!  incidence to the subgraph spanned by the ordered vertices.
!  Among all the vertices of maximal incidence, a vertex of
!  maximal degree is chosen. This subroutine determines the
!  inverse of the incidence degree ordering, that is, an array
!  listp such that listp(list(k)) = k for k = 1,2,...,n.

!  The subroutine statement is

!    subroutine idog(n, nghbrp, npntrp, nghbrs, npntrs, listp,
!                    maxclq, maxid, iwa1, iwa2, iwa3)

!  where

!    n is a positive integer input variable set to the number of vertices of G.

!    nghbrp is an integer input array which contains the
!      predecessor adjacency lists for the graph G.

!    npntrp is an integer input array of length n + 1 which specifies the
!      locations of the predecessor adjacency lists in nghbrp.
!      The vertices preceding and adjacent to vertex j are

!            nghbrp(k), k = npntrp(j),...,npntrp(j+1)-1.

!      Note that npntrp(n+1)-1 is then the number of vertices
!      plus edges of the graph G.

!    nghbrs is an integer input array which contains the
!      successor adjacency lists for the graph G.

!    npntrs is an integer input array of length n + 1 which specifies the
!      locations of the successor adjacency lists in nghbrs.
!      The vertices succeeding and adjacent to vertex j are

!            nghbrs(k), k = npntrs(j),...,npntrs(j+1)-1.

!      Note that npntrs(n+1)-1 is then the number of vertices
!      plus edges of the graph G.

!    listp is an integer output array of length n which specifies
!      the inverse of the incidence degree ordering of the
!      vertices.  Vertex j is in position listp(j) of this ordering.

!    maxclq is an integer output variable set to the size
!      of the largest clique found during the ordering.

!    maxid is an integer output variable set to the maximum
!      incidence degree found during the ordering.

!    iwa1, iwa2, and iwa3 are integer work arrays of length n.

!  Subprograms called

!    MINPACK-supplied ... numsrt

!    FORTRAN-supplied ... max

!  Argonne National Laboratory. MINPACK Project. December 1984.
!  Thomas F. Coleman, Burton S. Garbow, Jorge J. More'

!  **********

INTEGER :: i, j, k, maxinc, maxdeg, maxlst, ncomp, numdeg, numinc, numord,  &
           iwa1(0:n-1), iwa2(n), iwa3(n)

!     Initialization block.

DO  j = 1, n
  listp(j) = (npntrp(j+1) - npntrp(j) - 1) + (npntrs(j+1) - npntrs(j) - 1)
END DO
maxlst = (npntrp(n+1) + npntrs(n+1))/n

!     Sort the degree sequence.

CALL numsrt(n, n-1, listp, 1, iwa1, iwa2, iwa3)

!     Create a doubly-linked list to access the incidences of the
!     vertices. The pointers for the linked list are as follows.

!     Each un-ordered vertex i is in a list (the incidence list)
!     of vertices with the same incidence.

!     iwa1(numinc) is the first vertex in the numinc list unless
!     iwa1(numinc) = 0.  In this case there are no vertices in the numinc list.

!     iwa2(i) is the vertex before i in the incidence list unless iwa2(i) = 0.
!     In this case i is the first vertex in this incidence list.

!     iwa3(i) is the vertex after i in the incidence list unless iwa3(i) = 0.
!     In this case i is the last vertex in this incidence list.

!     If i is an un-ordered vertex, then -listp(i) is the incidence of i to the
!     graph induced by the ordered vertices.  If j is an ordered vertex, then
!     listp(j) is the incidence degree order of vertex j.

maxinc = 0
DO  j = 1, n
  i = iwa1(j-1)
  iwa1(j-1) = 0
  iwa2(i) = 0
  iwa3(i) = iwa1(0)
  IF (iwa1(0) > 0) iwa2(iwa1(0)) = i
  iwa1(0) = i
  listp(j) = 0
END DO
maxclq = 0
maxid = 0
numord = 1

!     Beginning of iteration loop.

!        Choose a vertex j of maximal degree among the
!        vertices of maximal incidence maxinc.

30 DO
  k = iwa1(maxinc)
  IF (k > 0) EXIT
  maxinc = maxinc - 1
END DO

maxdeg = -1
DO  i = 1, maxlst
  numdeg = (npntrp(k+1) - npntrp(k) - 1) + (npntrs(k+1) - npntrs(k) - 1)
  IF (numdeg > maxdeg) THEN
    maxdeg = numdeg
    j = k
  END IF
  k = iwa3(k)
  IF (k <= 0) EXIT
END DO

listp(j) = numord
maxid = MAX(maxid, maxinc)

!        Update the size of the largest clique found during the ordering.

IF (maxinc == 0) ncomp = 0
ncomp = ncomp + 1
IF (maxinc+1 == ncomp) maxclq = MAX(maxclq, ncomp)

!        Termination test.

numord = numord + 1
IF (numord > n) GO TO 100

!        Delete vertex j from the maxinc list.

IF (iwa2(j) == 0) THEN
  iwa1(maxinc) = iwa3(j)
ELSE
  iwa3(iwa2(j)) = iwa3(j)
END IF
IF (iwa3(j) > 0) iwa2(iwa3(j)) = iwa2(j)

!        Determine all the neighbors of vertex j which precede j
!        in the subgraph spanned by the un-ordered vertices.

DO  k = npntrp(j), npntrp(j+1)-1
  i = nghbrp(k)
  
!           Update the pointers to the current incidence lists.
  
  numinc = -listp(i)
  IF (numinc >= 0) THEN
    listp(i) = listp(i) - 1
    maxinc = MAX(maxinc,-listp(i))
    
!              Delete vertex i from the numinc list.
    
    IF (iwa2(i) == 0) THEN
      iwa1(numinc) = iwa3(i)
    ELSE
      iwa3(iwa2(i)) = iwa3(i)
    END IF
    IF (iwa3(i) > 0) iwa2(iwa3(i)) = iwa2(i)
    
!              Add vertex i to the numinc+1 list.
    
    iwa2(i) = 0
    iwa3(i) = iwa1(numinc+1)
    IF (iwa1(numinc+1) > 0) iwa2(iwa1(numinc+1)) = i
    iwa1(numinc+1) = i
  END IF
END DO

!        Determine all the neighbors of vertex j which succeed j
!        in the subgraph spanned by the un-ordered vertices.

DO  k = npntrs(j), npntrs(j+1)-1
  i = nghbrs(k)
  
!           Update the pointers to the current incidence lists.
  
  numinc = -listp(i)
  IF (numinc >= 0) THEN
    listp(i) = listp(i) - 1
    maxinc = MAX(maxinc,-listp(i))
    
!              Delete vertex i from the numinc list.
    
    IF (iwa2(i) == 0) THEN
      iwa1(numinc) = iwa3(i)
    ELSE
      iwa3(iwa2(i)) = iwa3(i)
    END IF
    IF (iwa3(i) > 0) iwa2(iwa3(i)) = iwa2(i)
    
!              Add vertex i to the numinc+1 list.
    
    iwa2(i) = 0
    iwa3(i) = iwa1(numinc+1)
    IF (iwa1(numinc+1) > 0) iwa2(iwa1(numinc+1)) = i
    iwa1(numinc+1) = i
  END IF
END DO

!        End of iteration loop.

GO TO 30
100 RETURN

!     Last card of subroutine idog.

END SUBROUTINE idog



SUBROUTINE numsrt(n, nmax, num, mode, INDEX, last, next)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-30  Time: 12:49:10

INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: nmax
INTEGER, INTENT(IN)   :: num(:)
INTEGER, INTENT(IN)   :: mode
INTEGER, INTENT(OUT)  :: INDEX(:)
INTEGER, INTENT(OUT)  :: last(0:)
INTEGER, INTENT(OUT)  :: next(:)

!  **********.

!  subroutine numsrt

!  Given a sequence of integers, this subroutine groups together those indices
!  with the same sequence value and, optionally, sorts the sequence into either
!  ascending or descending order.

!  The sequence of integers is defined by the array num, and it is assumed
!  that the integers are each from the set 0,1,...,nmax.  On output the indices
!  k such that num(k) = l for any l = 0,1,...,nmax can be obtained from the
!  arrays last and next as follows.

!        k = last(l)
!        while (k .ne. 0) k = next(k)

!  Optionally, the subroutine produces an array index so that
!  the sequence num(index(i)), i = 1,2,...,n is sorted.

!  The subroutine statement is

!    subroutine numsrt(n, nmax, num, mode, index, last, next)

!  where

!    n is a positive integer input variable.

!    nmax is a positive integer input variable.

!    num is an input array of length n which contains the sequence of
!      integers to be grouped and sorted.  It is assumed that the integers
!      are each from the set 0,1,...,nmax.

!    mode is an integer input variable.  The sequence num is sorted in
!      ascending order if mode is positive and in descending order if mode is
!      negative.  If mode is 0, no sorting is done.

!    index is an integer output array of length n set so that the sequence

!            num(index(i)), i = 1,2,...,n

!      is sorted according to the setting of mode.  If mode is 0,
!      index is not referenced.

!    last is an integer output array of length nmax + 1.  The index of num for
!      the last occurrence of l is last(l) for any l = 0,1,...,nmax unless
!      last(l) = 0.  In this case l does not appear in num.

!    next is an integer output array of length n.  If num(k) = l, then the
!      index of num for the previous occurrence of l is next(k) for any
!      l = 0,1,...,nmax unless next(k) = 0.
!      In this case there is no previous occurrence of l in num.

!  Argonne National Laboratory. MINPACK Project. July 1983.
!  Thomas F. Coleman, Burton S. Garbow, Jorge J. More'

!  **********

INTEGER :: i, j, jinc, jl, ju, k, l

!     Determine the arrays next and last.

last(0:nmax) = 0
DO  k = 1, n
  l = num(k)
  next(k) = last(l)
  last(l) = k
END DO
IF (mode == 0) RETURN

!     Store the pointers to the sorted array in index.

i = 1
IF (mode > 0) THEN
  jl = 0
  ju = nmax
  jinc = 1
ELSE
  jl = nmax
  ju = 0
  jinc = -1
END IF
DO  j = jl, ju, jinc
  k = last(j)
  DO
    IF (k == 0) EXIT
    INDEX(i) = k
    i = i + 1
    k = next(k)
  END DO
END DO
RETURN

!     Last card of subroutine numsrt.

END SUBROUTINE numsrt



SUBROUTINE sdpt(n, nghbrp, npntrp, nghbrs, npntrs, ngrp, maxgrp)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-30  Time: 12:49:17

INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: nghbrp(:)
INTEGER, INTENT(IN)   :: npntrp(:)    ! npntrp(n+1)
INTEGER, INTENT(IN)   :: nghbrs(:)
INTEGER, INTENT(IN)   :: npntrs(:)    ! npntrs(n+1)
INTEGER, INTENT(OUT)  :: ngrp(:)
INTEGER, INTENT(OUT)  :: maxgrp

!  **********

!  subroutine sdpt

!  Given a loopless graph G = (V,E), this subroutine determines
!  a symmetric coloring of G by the Powell-Toint direct method.

!  The Powell-Toint method assigns the k-th color by examining the un-colored
!  vertices U(k) in order of non-increasing degree and assigning color k to
!  vertex v if there are no paths of length 1 or 2 (in the graph induced
!  by U(k)) between v and some k-colored vertex.

!  The subroutine statement is

!    subroutine sdpt(n, nghbrp, npntrp, nghbrs, npntrs, ngrp, maxgrp,
!                    iwa1, iwa2)

!  where

!    n is a positive integer input variable set to the number of vertices of G.

!    nghbrp is an integer input array which contains the
!      predecessor adjacency lists for the graph G.

!    npntrp is an integer input array of length n + 1 which specifies the
!      locations of the predecessor adjacency lists in nghbrp.
!      The vertices preceding and adjacent to vertex j are

!            nghbrp(k), k = npntrp(j),...,npntrp(j+1)-1.

!      Note that npntrp(n+1)-1 is then the number of vertices
!      plus edges of the graph G.

!    nghbrs is an integer input array which contains the
!      successor adjacency lists for the graph G.

!    npntrs is an integer input array of length n + 1 which specifies the
!      locations of the successor adjacency lists in nghbrs.
!      The vertices succeeding and adjacent to vertex j are

!            nghbrs(k), k = npntrs(j),...,npntrs(j+1)-1.

!      Note that npntrs(n+1)-1 is then the number of vertices
!      plus edges of the graph G.

!    ngrp is an integer output array of length n which specifies the
!      symmetric coloring of G. Vertex j is colored with color ngrp(j).

!    maxgrp is an integer output variable which specifies the
!      number of colors in the symmetric coloring of G.

!    iwa1 and iwa2 are integer work arrays of length n.

!  Subprograms called

!    FORTRAN-supplied ... max

!  Argonne National Laboratory. MINPACK Project. December 1984.
!  Thomas F. Coleman, Burton S. Garbow, Jorge J. More'

!  **********

INTEGER :: j, jp, k, kp, l, maxdeg, numdeg, numv, iwa1(0:n-1), iwa2(n)

!     Initialization block. Numv is the current number of un-colored
!     vertices, maxdeg is the maximum induced degree of these
!     vertices, and maxgrp is the current group number (color).

numv = n
maxdeg = 0
DO  j = 1, n
  ngrp(j) = (npntrp(j) - npntrp(j+1) + 1) + (npntrs(j) - npntrs(j+1) + 1)
  maxdeg = MAX(maxdeg,-ngrp(j))
  iwa2(j) = -j
END DO
maxgrp = 0

!     Beginning of iteration loop.

!        Sort the list of un-colored vertices so that their
!        induced degrees are in non-decreasing order.

20 iwa1(0:maxdeg) = 0
DO  l = 1, numv
  numdeg = -ngrp(-iwa2(l))
  iwa1(numdeg) = iwa1(numdeg) + 1
END DO
k = 1
DO  numdeg = maxdeg, 0, -1
  l = iwa1(numdeg)
  iwa1(numdeg) = k
  k = k + l
END DO
k = 1

60 j = iwa2(k)
IF (j > 0) THEN
  k = iwa1(-ngrp(j))
ELSE
  numdeg = -ngrp(-j)
  l = iwa1(numdeg)
  iwa2(k) = iwa2(l)
  iwa2(l) = -j
  iwa1(numdeg) = iwa1(numdeg) + 1
END IF
IF (k <= numv) GO TO 60
maxgrp = maxgrp + 1

!        Determine the vertices in group maxgrp.

loop150:   &
DO  l = 1, numv
  j = iwa2(l)
  
!          Examine each vertex k preceding vertex j and all the neighbors of
!          vertex k to determine if vertex j can be considered for group maxgrp.
  
  DO  jp = npntrp(j), npntrp(j+1)-1
    k = nghbrp(jp)
    IF (ngrp(k) == maxgrp) CYCLE loop150
    IF (ngrp(k) <= 0) THEN
      DO  kp = npntrp(k), npntrp(k+1)-1
        IF (ngrp(nghbrp(kp)) == maxgrp) CYCLE loop150
      END DO
      DO  kp = npntrs(k), npntrs(k+1)-1
        IF (ngrp(nghbrs(kp)) == maxgrp) CYCLE loop150
      END DO
    END IF
  END DO
  
!          Examine each vertex k succeeding vertex j and all the neighbors of
!          vertex k to determine if vertex j can be added to group maxgrp.
  
  DO  jp = npntrs(j), npntrs(j+1)-1
    k = nghbrs(jp)
    IF (ngrp(k) == maxgrp) CYCLE loop150
    IF (ngrp(k) <= 0) THEN
      DO  kp = npntrp(k), npntrp(k+1)-1
        IF (ngrp(nghbrp(kp)) == maxgrp) CYCLE loop150
      END DO
      DO  kp = npntrs(k), npntrs(k+1)-1
        IF (ngrp(nghbrs(kp)) == maxgrp) CYCLE loop150
      END DO
    END IF
  END DO
  
!           Add vertex j to group maxgrp and remove vertex j
!           from the list of un-colored vertices.
  
  ngrp(j) = maxgrp
  iwa2(l) = 0
  
!           Update the degrees of the neighbors of vertex j.
  
  DO  jp = npntrp(j), npntrp(j+1)-1
    k = nghbrp(jp)
    IF (ngrp(k) < 0) ngrp(k) = ngrp(k) + 1
  END DO
  DO  jp = npntrs(j), npntrs(j+1)-1
    k = nghbrs(jp)
    IF (ngrp(k) < 0) ngrp(k) = ngrp(k) + 1
  END DO

END DO loop150

!        Compress the updated list of un-colored vertices.
!        Reset numv and recompute maxdeg.

k = 0
maxdeg = 0
DO  l = 1, numv
  IF (iwa2(l) /= 0) THEN
    k = k + 1
    iwa2(k) = -iwa2(l)
    maxdeg = MAX(maxdeg,-ngrp(iwa2(l)))
  END IF
END DO
numv = k

!        End of iteration loop.

IF (numv > 0) GO TO 20
RETURN

!     Last card of subroutine sdpt.

END SUBROUTINE sdpt



SUBROUTINE seq(n, indrow, jpntr, indcol, ipntr, list, ngrp, maxgrp)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-30  Time: 12:49:22

INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: indrow(:)
INTEGER, INTENT(IN)   :: jpntr(n+1)
INTEGER, INTENT(IN)   :: indcol(:)
INTEGER, INTENT(IN)   :: ipntr(:)
INTEGER, INTENT(IN)   :: list(:)
INTEGER, INTENT(OUT)  :: ngrp(:)
INTEGER, INTENT(OUT)  :: maxgrp

!  **********

!  subroutine seq

!  Given the sparsity pattern of an m by n matrix A, this
!  subroutine determines a consistent partition of the
!  columns of A by a sequential algorithm.

!  A consistent partition is defined in terms of the loopless
!  graph G with vertices a(j), j = 1,2,...,n where a(j) is the
!  j-th column of A and with edge (a(i),a(j)) if and only if
!  columns i and j have a non-zero in the same row position.

!  A partition of the columns of A into groups is consistent
!  if the columns in any group are not adjacent in the graph G.
!  In graph-theory terminology, a consistent partition of the
!  columns of A corresponds to a coloring of the graph G.

!  The subroutine examines the columns in the order specified
!  by the array list, and assigns the current column to the
!  group with the smallest possible number.

!  Note that the value of m is not needed by seq and is
!  therefore not present in the subroutine statement.

!  The subroutine statement is

!    subroutine seq(n, indrow, jpntr, indcol, ipntr, list, ngrp, maxgrp, iwa)

!  where

!    n is a positive integer input variable set to the number of columns of A.

!    indrow is an integer input array which contains the row
!      indices for the non-zeroes in the matrix A.

!    jpntr is an integer input array of length n + 1 which
!      specifies the locations of the row indices in indrow.
!      The row indices for column j are

!            indrow(k), k = jpntr(j),...,jpntr(j+1)-1.

!      Note that jpntr(n+1)-1 is then the number of non-zero
!      elements of the matrix A.

!    indcol is an integer input array which contains the
!      column indices for the non-zeroes in the matrix A.

!    ipntr is an integer input array of length m + 1 which
!      specifies the locations of the column indices in indcol.
!      The column indices for row i are

!            indcol(k), k = ipntr(i),...,ipntr(i+1)-1.

!      Note that ipntr(m+1)-1 is then the number of non-zero
!      elements of the matrix A.

!    list is an integer input array of length n which specifies
!      the order to be used by the sequential algorithm.
!      The j-th column in this order is list(j).

!    ngrp is an integer output array of length n which specifies the
!      partition of the columns of A.  Column jcol belongs to group ngrp(jcol).

!    maxgrp is an integer output variable which specifies the
!      number of groups in the partition of the columns of A.

!    iwa is an integer work array of length n.

!  Argonne National Laboratory. MINPACK Project. July 1983.
!  Thomas F. Coleman, Burton S. Garbow, Jorge J. More'

!  **********

INTEGER :: ic, ip, ir, j, jcol, jp, iwa(n)

!     Initialization block.

maxgrp = 0
DO  jp = 1, n
  ngrp(jp) = n
  iwa(jp) = 0
END DO

!     Beginning of iteration loop.

DO  j = 1, n
  jcol = list(j)
  
!        Find all columns adjacent to column jcol.
  
!        Determine all positions (ir,jcol) which correspond
!        to non-zeroes in the matrix.
  
  DO  jp = jpntr(jcol), jpntr(jcol+1)-1
    ir = indrow(jp)
    
!           For each row ir, determine all positions (ir,ic)
!           which correspond to non-zeroes in the matrix.
    
    DO  ip = ipntr(ir), ipntr(ir+1)-1
      ic = indcol(ip)
      
!              Array iwa marks the group numbers of the
!              columns which are adjacent to column jcol.
      
      iwa(ngrp(ic)) = j
    END DO
  END DO
  
!        Assign the smallest un-marked group number to jcol.
  
  DO  jp = 1, maxgrp
    IF (iwa(jp) /= j) GO TO 50
  END DO
  maxgrp = maxgrp + 1
  50 ngrp(jcol) = jp
END DO

!        End of iteration loop.

RETURN

!     Last card of subroutine seq.

END SUBROUTINE seq



SUBROUTINE setr(m, n, indrow, jpntr, indcol, ipntr)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-30  Time: 12:49:30

INTEGER, INTENT(IN)   :: m
INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: indrow(:)
INTEGER, INTENT(IN)   :: jpntr(:)
INTEGER, INTENT(OUT)  :: indcol(:)
INTEGER, INTENT(OUT)  :: ipntr(:)

!  **********

!  subroutine setr

!  Given a column-oriented definition of the sparsity pattern
!  of an m by n matrix A, this subroutine determines a
!  row-oriented definition of the sparsity pattern of A.

!  On input the column-oriented definition is specified by
!  the arrays indrow and jpntr. On output the row-oriented
!  definition is specified by the arrays indcol and ipntr.

!  The subroutine statement is

!    subroutine setr(m, n, indrow, jpntr, indcol, ipntr, iwa)

!  where

!    m is a positive integer input variable set to the number of rows of A.

!    n is a positive integer input variable set to the number of columns of A.

!    indrow is an integer input array which contains the row
!      indices for the non-zeroes in the matrix A.

!    jpntr is an integer input array of length n + 1 which
!      specifies the locations of the row indices in indrow.
!      The row indices for column j are

!            indrow(k), k = jpntr(j),...,jpntr(j+1)-1.

!      Note that jpntr(n+1)-1 is then the number of non-zero
!      elements of the matrix A.

!    indcol is an integer output array which contains the
!      column indices for the non-zeroes in the matrix A.

!    ipntr is an integer output array of length m + 1 which
!      specifies the locations of the column indices in indcol.
!      The column indices for row i are

!            indcol(k), k = ipntr(i),...,ipntr(i+1)-1.

!      Note that ipntr(1) is set to 1 and that ipntr(m+1)-1 is
!      then the number of non-zero elements of the matrix A.

!    iwa is an integer work array of length m.

!  Argonne National Laboratory. MINPACK Project. July 1983.
!  Thomas F. Coleman, Burton S. Garbow, Jorge J. More'

!  **********

INTEGER :: ir, jcol, jp, iwa(n)

!     Store in array iwa the counts of non-zeroes in the rows.

iwa(1:m) = 0
DO  jp = 1, jpntr(n+1)-1
  iwa(indrow(jp)) = iwa(indrow(jp)) + 1
END DO

!     Set pointers to the start of the rows in indcol.

ipntr(1) = 1
DO  ir = 1, m
  ipntr(ir+1) = ipntr(ir) + iwa(ir)
  iwa(ir) = ipntr(ir)
END DO

!     Fill indcol.

DO  jcol = 1, n
  DO  jp = jpntr(jcol), jpntr(jcol+1)-1
    ir = indrow(jp)
    indcol(iwa(ir)) = jcol
    iwa(ir) = iwa(ir) + 1
  END DO
END DO
RETURN

!     Last card of subroutine setr.

END SUBROUTINE setr



SUBROUTINE slo(n, indrow, jpntr, indcol, ipntr, ndeg, list, maxclq)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-30  Time: 12:49:41

INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: indrow(:)
INTEGER, INTENT(IN)   :: jpntr(:)
INTEGER, INTENT(IN)   :: indcol(:)
INTEGER, INTENT(IN)   :: ipntr(:)
INTEGER, INTENT(IN)   :: ndeg(:)
INTEGER, INTENT(OUT)  :: list(:)
INTEGER, INTENT(OUT)  :: maxclq

!  **********

!  subroutine slo

!  Given the sparsity pattern of an m by n matrix A, this subroutine
!  determines the smallest-last ordering of the columns of A.

!  The smallest-last ordering is defined for the loopless
!  graph G with vertices a(j), j = 1,2,...,n where a(j) is the
!  j-th column of A and with edge (a(i),a(j)) if and only if
!  columns i and j have a non-zero in the same row position.

!  The smallest-last ordering is determined recursively by
!  letting list(k), k = n,...,1 be a column with least degree
!  in the subgraph spanned by the un-ordered columns.

!  Note that the value of m is not needed by slo and is
!  therefore not present in the subroutine statement.

!  The subroutine statement is

!    subroutine slo(n, indrow, jpntr, indcol, ipntr, ndeg, list,
!                   maxclq, iwa1, iwa2, iwa3, iwa4)

!  where

!    n is a positive integer input variable set to the number of columns of A.

!    indrow is an integer input array which contains the row
!      indices for the non-zeroes in the matrix A.

!    jpntr is an integer input array of length n + 1 which
!      specifies the locations of the row indices in indrow.
!      The row indices for column j are

!            indrow(k), k = jpntr(j),...,jpntr(j+1)-1.

!      Note that jpntr(n+1)-1 is then the number of non-zero
!      elements of the matrix A.

!    indcol is an integer input array which contains the
!      column indices for the non-zeroes in the matrix A.

!    ipntr is an integer input array of length m + 1 which
!      specifies the locations of the column indices in indcol.
!      The column indices for row i are

!            indcol(k), k = ipntr(i),...,ipntr(i+1)-1.

!      Note that ipntr(m+1)-1 is then the number of non-zero
!      elements of the matrix A.

!    ndeg is an integer input array of length n which specifies
!      the degree sequence. The degree of the j-th column
!      of A is ndeg(j).

!    list is an integer output array of length n which specifies
!      the smallest-last ordering of the columns of A.  The j-th
!      column in this order is list(j).

!    maxclq is an integer output variable set to the size
!      of the largest clique found during the ordering.

!    iwa1,iwa2,iwa3, and iwa4 are integer work arrays of length n.

!  Subprograms called

!    FORTRAN-supplied ... min

!  Argonne National Laboratory. MINPACK Project. August 1984.
!  Thomas F. Coleman, Burton S. Garbow, Jorge J. More'

!  **********

INTEGER :: ic, ip, ir, jcol, jp, mindeg, numdeg, numord,  &
           iwa1(0:n-1), iwa2(n), iwa3(n), iwa4(n)

!     Initialization block.

mindeg = n
DO  jp = 1, n
  iwa1(jp-1) = 0
  iwa4(jp) = n
  list(jp) = ndeg(jp)
  mindeg = MIN(mindeg, ndeg(jp))
END DO

!     Create a doubly-linked list to access the degrees of the
!     columns. The pointers for the linked list are as follows.

!     Each un-ordered column ic is in a list (the degree list)
!     of columns with the same degree.

!     iwa1(numdeg) is the first column in the numdeg list unless
!     iwa1(numdeg) = 0.  In this case there are no columns in the numdeg list.

!     iwa2(ic) is the column before ic in the degree list unless iwa2(ic) = 0.
!     In this case ic is the first column in this degree list.

!     iwa3(ic) is the column after ic in the degree list unless iwa3(ic) = 0.
!     In this case ic is the last column in this degree list.

!     If ic is an un-ordered column, then list(ic) is the degree of ic in the
!     graph induced by the un-ordered columns.  If jcol is an ordered column,
!     then list(jcol) is the smallest-last order of column jcol.

DO  jp = 1, n
  numdeg = ndeg(jp)
  iwa2(jp) = 0
  iwa3(jp) = iwa1(numdeg)
  IF (iwa1(numdeg) > 0) iwa2(iwa1(numdeg)) = jp
  iwa1(numdeg) = jp
END DO
maxclq = 0
numord = n

!     Beginning of iteration loop.

!        Choose a column jcol of minimal degree mindeg.

30 DO
  jcol = iwa1(mindeg)
  IF (jcol > 0) EXIT
  mindeg = mindeg + 1
END DO

list(jcol) = numord

!        Mark the size of the largest clique found during the ordering.

IF (mindeg+1 == numord .AND. maxclq == 0) maxclq = numord

!        Termination test.

numord = numord - 1
IF (numord == 0) GO TO 80

!        Delete column jcol from the mindeg list.

iwa1(mindeg) = iwa3(jcol)
IF (iwa3(jcol) > 0) iwa2(iwa3(jcol)) = 0

!        Find all columns adjacent to column jcol.

iwa4(jcol) = 0

!        Determine all positions (ir,jcol) which correspond to non-zeroes
!        in the matrix.

DO  jp = jpntr(jcol), jpntr(jcol+1)-1
  ir = indrow(jp)
  
!           For each row ir, determine all positions (ir,ic)
!           which correspond to non-zeroes in the matrix.
  
  DO  ip = ipntr(ir), ipntr(ir+1)-1
    ic = indcol(ip)
    
!              Array iwa4 marks columns which are adjacent to column jcol.
    
    IF (iwa4(ic) > numord) THEN
      iwa4(ic) = numord
      
!                 Update the pointers to the current degree lists.
      
      numdeg = list(ic)
      list(ic) = list(ic) - 1
      mindeg = MIN(mindeg,list(ic))
      
!                 Delete column ic from the numdeg list.
      
      IF (iwa2(ic) == 0) THEN
        iwa1(numdeg) = iwa3(ic)
      ELSE
        iwa3(iwa2(ic)) = iwa3(ic)
      END IF
      IF (iwa3(ic) > 0) iwa2(iwa3(ic)) = iwa2(ic)
      
!                 Add column ic to the numdeg-1 list.
      
      iwa2(ic) = 0
      iwa3(ic) = iwa1(numdeg-1)
      IF (iwa1(numdeg-1) > 0) iwa2(iwa1(numdeg-1)) = ic
      iwa1(numdeg-1) = ic
    END IF
  END DO
END DO

!        End of iteration loop.

GO TO 30

!     Invert the array list.

80 DO  jcol = 1, n
  iwa2(list(jcol)) = jcol
END DO
list(1:n) = iwa2(1:n)
RETURN

!     Last card of subroutine slo.

END SUBROUTINE slo



SUBROUTINE slog(n, nghbrp, npntrp, nghbrs, npntrs, listp, maxclq, maxvd)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-30  Time: 12:49:46

INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: nghbrp(:)
INTEGER, INTENT(IN)   :: npntrp(:)
INTEGER, INTENT(IN)   :: nghbrs(:)
INTEGER, INTENT(IN)   :: npntrs(:)
INTEGER, INTENT(OUT)  :: listp(:)
INTEGER, INTENT(OUT)  :: maxclq
INTEGER, INTENT(OUT)  :: maxvd

!  **********

!  subroutine slog

!  Given a loopless graph G = (V,E), this subroutine determines
!  the smallest-last ordering of the vertices of G.

!  The smallest-last ordering is determined recursively by
!  letting list(k), k = n,...,1 be a vertex with least degree
!  in the subgraph spanned by the un-ordered vertices.
!  This subroutine determines the inverse of the smallest-last
!  ordering, that is, an array listp such that listp(list(k)) = k
!  for k = 1,2,...,n.

!  The subroutine statement is

!    subroutine slog(n, nghbrp, npntrp, nghbrs, npntrs, listp,
!                    maxclq, maxvd, iwa1, iwa2, iwa3)

!  where

!    n is a positive integer input variable set to the number
!      of vertices of G.

!    nghbrp is an integer input array which contains the
!      predecessor adjacency lists for the graph G.

!    npntrp is an integer input array of length n + 1 which specifies
!      the locations of the predecessor adjacency lists in nghbrp.
!      The vertices preceding and adjacent to vertex j are

!            nghbrp(k), k = npntrp(j),...,npntrp(j+1)-1.

!      Note that npntrp(n+1)-1 is then the number of vertices
!      plus edges of the graph G.

!    nghbrs is an integer input array which contains the
!      successor adjacency lists for the graph G.

!    npntrs is an integer input array of length n + 1 which specifies
!      the locations of the successor adjacency lists in nghbrs.
!      The vertices succeeding and adjacent to vertex j are

!            nghbrs(k), k = npntrs(j),...,npntrs(j+1)-1.

!      Note that npntrs(n+1)-1 is then the number of vertices
!      plus edges of the graph G.

!    listp is an integer output array of length n which specifies
!      the inverse of the smallest-last ordering of the vertices.
!      Vertex j is in position listp(j) of this ordering.

!    maxclq is an integer output variable set to the size
!      of the largest clique found during the ordering.

!    maxvd is an integer output variable set to the maximum
!      vertex degree found during the ordering.

!    iwa1,iwa2, and iwa3 are integer work arrays of length n.

!  Subprograms called

!    FORTRAN-supplied ... max,min

!  Argonne National Laboratory. MINPACK Project. December 1984.
!  Thomas F. Coleman, Burton S. Garbow, Jorge J. More'

!  **********

INTEGER :: i, j, k, mindeg, numdeg, numord, iwa1(0:n-1), iwa2(n), iwa3(n)

!     Initialization block.

mindeg = n
DO  j = 1, n
  iwa1(j-1) = 0
  listp(j) = (npntrp(j) - npntrp(j+1) + 1) + (npntrs(j) - npntrs(j+1) + 1)
  mindeg = MIN(mindeg,-listp(j))
END DO

!     Create a doubly-linked list to access the degrees of the
!     vertices.  The pointers for the linked list are as follows.

!     Each un-ordered vertex i is in a list (the degree list)
!     of vertices with the same degree.

!     iwa1(numdeg) is the first vertex in the numdeg list unless
!     iwa1(numdeg) = 0.  In this case there are no vertices in the numdeg list.

!     iwa2(i) is the vertex before i in the degree list unless iwa2(i) = 0.
!     In this case i is the first vertex in this degree list.

!     iwa3(i) is the vertex after i in the degree list unless iwa3(i) = 0.
!     In this case i is the last vertex in this degree list.

!     If i is an un-ordered vertex, then -listp(i) is the degree of i in the
!     graph induced by the un-ordered vertices.  If j is an ordered vertex,
!     then listp(j) is the smallest-last order of vertex j.

DO  j = 1, n
  numdeg = -listp(j)
  iwa2(j) = 0
  iwa3(j) = iwa1(numdeg)
  IF (iwa1(numdeg) > 0) iwa2(iwa1(numdeg)) = j
  iwa1(numdeg) = j
END DO
maxclq = 0
maxvd = 0
numord = n

!     Beginning of iteration loop.

!        Choose a vertex j of minimal degree mindeg.

30 DO
  j = iwa1(mindeg)
  IF (j > 0) EXIT
  mindeg = mindeg + 1
END DO

listp(j) = numord
maxvd = MAX(maxvd, mindeg)

!        Mark the size of the largest clique
!        found during the ordering.

IF (mindeg+1 == numord .AND. maxclq == 0) maxclq = numord

!        Termination test.

numord = numord - 1
IF (numord == 0) GO TO 80

!        Delete vertex j from the mindeg list.

iwa1(mindeg) = iwa3(j)
IF (iwa3(j) > 0) iwa2(iwa3(j)) = 0

!        Determine all the neighbors of vertex j which precede j
!        in the subgraph spanned by the un-ordered vertices.

DO  k = npntrp(j), npntrp(j+1)-1
  i = nghbrp(k)
  
!           Update the pointers to the current degree lists.
  
  numdeg = -listp(i)
  IF (numdeg >= 0) THEN
    listp(i) = listp(i) + 1
    mindeg = MIN(mindeg,-listp(i))
    
!              Delete vertex i from the numdeg list.
    
    IF (iwa2(i) == 0) THEN
      iwa1(numdeg) = iwa3(i)
    ELSE
      iwa3(iwa2(i)) = iwa3(i)
    END IF
    IF (iwa3(i) > 0) iwa2(iwa3(i)) = iwa2(i)
    
!              Add vertex i to the numdeg-1 list.
    
    iwa2(i) = 0
    iwa3(i) = iwa1(numdeg-1)
    IF (iwa1(numdeg-1) > 0) iwa2(iwa1(numdeg-1)) = i
    iwa1(numdeg-1) = i
  END IF
END DO

!        Determine all the neighbors of vertex j which succeed j
!        in the subgraph spanned by the un-ordered vertices.

DO  k = npntrs(j), npntrs(j+1)-1
  i = nghbrs(k)
  
!           Update the pointers to the current degree lists.
  
  numdeg = -listp(i)
  IF (numdeg >= 0) THEN
    listp(i) = listp(i) + 1
    mindeg = MIN(mindeg,-listp(i))
    
!              Delete vertex i from the numdeg list.
    
    IF (iwa2(i) == 0) THEN
      iwa1(numdeg) = iwa3(i)
    ELSE
      iwa3(iwa2(i)) = iwa3(i)
    END IF
    IF (iwa3(i) > 0) iwa2(iwa3(i)) = iwa2(i)
    
!              Add vertex i to the numdeg-1 list.
    
    iwa2(i) = 0
    iwa3(i) = iwa1(numdeg-1)
    IF (iwa1(numdeg-1) > 0) iwa2(iwa1(numdeg-1)) = i
    iwa1(numdeg-1) = i
  END IF
END DO

!        End of iteration loop.

GO TO 30
80 RETURN

!     Last card of subroutine slog.

END SUBROUTINE slog



SUBROUTINE srtdat(n, nnz, indrow, indcol, jpntr)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-30  Time: 12:49:52

INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN)      :: nnz
INTEGER, INTENT(IN OUT)  :: indrow(:)
INTEGER, INTENT(IN OUT)  :: indcol(:)
INTEGER, INTENT(OUT)     :: jpntr(:)    ! jpntr(n+1)

!  **********

!  subroutine srtdat

!  Given the non-zero elements of an m by n matrix A in
!  arbitrary order as specified by their row and column
!  indices, this subroutine permutes these elements so
!  that their column indices are in non-decreasing order.

!  On input it is assumed that the elements are specified in

!        indrow(k),indcol(k), k = 1,...,nnz.

!  On output the elements are permuted so that indcol is
!  in non-decreasing order. In addition, the array jpntr
!  is set so that the row indices for column j are

!        indrow(k), k = jpntr(j),...,jpntr(j+1)-1.

!  Note that the value of m is not needed by srtdat and is
!  therefore not present in the subroutine statement.

!  The subroutine statement is

!    subroutine srtdat(n,nnz,indrow,indcol,jpntr,iwa)

!  where

!    n is a positive integer input variable set to the number of columns of A.

!    nnz is a positive integer input variable set to the number
!      of non-zero elements of A.

!    indrow is an integer array of length nnz. On input indrow
!      must contain the row indices of the non-zero elements of A.
!      On output indrow is permuted so that the corresponding
!      column indices of indcol are in non-decreasing order.

!    indcol is an integer array of length nnz. On input indcol
!      must contain the column indices of the non-zero elements
!      of A. On output indcol is permuted so that these indices
!      are in non-decreasing order.

!    jpntr is an integer output array of length n + 1 which
!      specifies the locations of the row indices in the output
!      indrow. The row indices for column j are

!            indrow(k), k = jpntr(j),...,jpntr(j+1)-1.

!      Note that jpntr(1) is set to 1 and that jpntr(n+1)-1 is then nnz.

!    iwa is an integer work array of length n.

!  Subprograms called

!    FORTRAN-supplied ... max

!  Argonne National Laboratory. MINPACK Project. July 1983.
!  Thomas F. Coleman, Burton S. Garbow, Jorge J. More'

!  **********

INTEGER :: i, j, k, l, iwa(n)

!     Store in array iwa the counts of non-zeroes in the columns.

iwa(1:n) = 0
DO  k = 1, nnz
  iwa(indcol(k)) = iwa(indcol(k)) + 1
END DO

!     Set pointers to the start of the columns in indrow.

jpntr(1) = 1
DO  j = 1, n
  jpntr(j+1) = jpntr(j) + iwa(j)
  iwa(j) = jpntr(j)
END DO
k = 1

!     Begin in-place sort.

40 j = indcol(k)
IF (k >= jpntr(j)) THEN
  
!           Current element is in position.  Now examine the next element
!           or the first un-sorted element in the j-th group.
  
  k = MAX(k+1,iwa(j))
ELSE
  
!           Current element is not in position.  Place element in position
!           and make the displaced element the current element.
  
  l = iwa(j)
  iwa(j) = iwa(j) + 1
  i = indrow(k)
  indrow(k) = indrow(l)
  indcol(k) = indcol(l)
  indrow(l) = i
  indcol(l) = j
END IF
IF (k <= nnz) GO TO 40
RETURN

!     Last card of subroutine srtdat.

END SUBROUTINE srtdat

END MODULE coloring
