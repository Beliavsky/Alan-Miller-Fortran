MODULE Common_SPGMA2

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

!     PARAMETERS
INTEGER, PARAMETER  :: npmax = 50000, nvsmax = 30, nmax = npmax*2

!     COMMON SCALARS
INTEGER, SAVE    :: np, totnvs

!     COMMON ARRAYS
REAL (dp), SAVE  :: edges(npmax*nvsmax*3), vert(npmax*nvsmax*2)
INTEGER, SAVE    :: nvs(npmax)


CONTAINS


SUBROUTINE projp(p, base, x, y, xproj, yproj, inform)

!  This subroutine complements proj subroutine projecting
!  a point z_i in R^2 onto its corresponding polygon P_i.

!  On Entry:

!  p     integer,
!        index of the point z_p to be projected,

!  base  integer,
!        index which indicates the positions where the edges of the polygon
!        P_p are stored inside the "edges" vector which is part of the
!        polygons description.

!  x     REAL (dp),
!        first coordinate of z_i, i.e., z_i^1,

!  y     REAL (dp),
!        second coordinate of z_i, i.e., z_i^2,

!  On Return

!  xproj REAL (dp),
!        the projection of x,

!  yproj REAL (dp),
!        the projection of y.

!  inform integer,
!        termination parameter:
!        0 = the projection was successfully done,
!        1 = some error occurs in the projection.

!     PARAMETERS

INTEGER, INTENT(IN)     :: p
INTEGER, INTENT(IN)     :: base
REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(IN)   :: y
REAL (dp), INTENT(OUT)  :: xproj
REAL (dp), INTENT(OUT)  :: yproj
INTEGER, INTENT(OUT)    :: inform

!     COMMON SCALARS
! INTEGER :: np, totnvs

!     COMMON ARRAYS
! REAL (dp)  :: edges(npmax*nvsmax*3), vert(npmax*nvsmax*2)
! INTEGER    :: nvs(npmax)

!     LOCAL SCALARS
REAL (dp)  :: a, b, c, dist, mindis, nearvx, nearvy, px, py,  &
              v1x, v1y, v2x, v2y, vx, vy
INTEGER    :: i
LOGICAL    :: fact

!     LOCAL ARRAYS
LOGICAL    :: satisf(nvsmax)

!     COMMON BLOCKS
! COMMON /polyg/ nvs, vert, edges, np, totnvs

!     TEST ALL EDGES OF THE POLYGON FOR SATISFIABILITY

fact = .true.
inform = 0
DO i = base + 1, base + nvs(p)
  a = edges(3*i-2)
  b = edges(3*i-1)
  c = edges(3*i)
  IF (a*x + b*y + c <= 0.0D0) THEN
    satisf(i-base) = .true.

  ELSE
    fact = .false.
    satisf(i-base) = .false.
  END IF
END DO

IF (fact) THEN
  xproj = x
  yproj = y
  RETURN

END IF

!     SEARCH FOR THE CLOSEST VERTEX

mindis = 1.0D+99
DO i = base + 1, base + nvs(p)
  vx = vert(2*i-1)
  vy = vert(2*i)
  dist = SQRT((vx-x)**2 + (vy-y)**2)
  IF (dist <= mindis) THEN
    mindis = dist
    nearvx = vx
    nearvy = vy
  END IF

END DO

xproj = nearvx
yproj = nearvy

!     PROJECT ONTO THE VIOLATED CONSTRAINTS

DO i = base + 1, base + nvs(p)
  IF (.NOT.satisf(i-base)) THEN
    a = edges(3*i-2)
    b = edges(3*i-1)
    c = edges(3*i)
    IF (a /= 0.0D0) THEN
      py = (y - b*(c/a+x)/a) / (b**2/a**2 + 1)
      px = -(c + b*py) / a

    ELSE
      IF (b /= 0.0D0) THEN
        px = (x - a*(c/b+y)/b) / (a**2/b**2+1)
        py = -(c + a*px) / b

      ELSE
        WRITE (*,FMT=*) 'ERROR IN PROBLEM DEFINITION (a=b=0 for a ',  &
                        'constraint of the type a x + b y + c = 0)'
        inform = 1
      END IF

    END IF

    v1x = vert(2*i-1)
    v1y = vert(2*i)
    IF (i /= base + nvs(p)) THEN
      v2x = vert(2*i+1)
      v2y = vert(2*i+2)

    ELSE
      v2x = vert(2*base+1)
      v2y = vert(2*base+2)
    END IF

    IF (MIN(v1x,v2x) <= px .AND. px <= MAX(v1x,v2x) .AND.  &
        MIN(v1y,v2y) <= py .AND. py <= MAX(v1y,v2y)) THEN
      dist = SQRT((px-x)**2 + (py-y)**2)
      IF (dist <= mindis) THEN
        mindis = dist
        xproj = px
        yproj = py
      END IF

    END IF

  END IF

END DO
RETURN
END SUBROUTINE projp



SUBROUTINE genpro(nx, ny, xstep, ystep, prob, nvsvmi, nvsvma, rmin, rmax)

!  This subroutine generates a location problem (see the Figure).
!  First, a regular grid with nx horizontal points and ny vertical points in
!  the positive orthant is considered.  The points of the grid start at the
!  origin with an horizontal distance of xstep and a vertical distance of
!  ystept.  This grid will be the space at which the cities (represented by
!  polygons) will be distributed.  Before building the cities, an area
!  (rectangle) of preservation where almost nothing can be done, is defined.
!  This area of preservation will receive, after the construction of the
!  cities, an hydraulic plant of energy generation (to supply the energy to
!  the cities).  Then, in the rest of the space, the cities are built.
!  At each point of the grid (out of the central region) a city (represented
!  by a polygon) will be built with probability prob.  The definition of the
!  polygon uses variables nvsvmi, nvsvma, rmin and rmax in a way described in
!  the genpol (generate polygon) subroutine.  To transmit the energy from the
!  plant to the cities, a tower inside each city and a tower inside the
!  central region must be built.  The objective of the problem is to
!  determine the location of this towers in order to minimize the sum of the
!  distances from each city tower to the central one.

!  On Entry:

!  nx    integer,
!        number of horizontal points in the grid,

!  ny    integer,
!        number of vertical points in the grid,

!  xstep REAL (dp),
!        horizontal distance between points of the grid,

!  ystep REAL (dp),
!        vertical distance between points of the grid,

!  prob  REAL (dp),
!        probability of defining a city at point of the grid
!        (0 <= prob <= 1),

!  nvsvmi integer,
!        parameter for the polygon generation (described in genpol subroutine),

!  nvsvma integer,
!        parameter for the polygon generation (described in genpol subroutine),

!  rmin  REAL (dp),
!        parameter for the polygon generation (described in genpol subroutine),

!  rmax  REAL (dp),
!        parameter for the polygon generation (described in genpol subroutine).

!  On output:

!  As described in the genpol subroutine, the output is saved in the polyg
!  common block.

!     PARAMETERS

INTEGER, INTENT(IN)        :: nx
INTEGER, INTENT(IN)        :: ny
REAL (dp), INTENT(IN)      :: xstep
REAL (dp), INTENT(IN)      :: ystep
REAL (dp), INTENT(IN OUT)  :: prob
INTEGER, INTENT(IN OUT)    :: nvsvmi
INTEGER, INTENT(IN OUT)    :: nvsvma
REAL (dp), INTENT(IN OUT)  :: rmin
REAL (dp), INTENT(IN OUT)  :: rmax

!     COMMON SCALARS
! INTEGER  :: np, totnvs

!     COMMON ARRAYS
! REAL (dp)  :: edges(npmax*nvsmax*3), vert(npmax*nvsmax*2)
! INTEGER    :: nvs(npmax)

!     LOCAL SCALARS
REAL (dp)  :: cx, cy, lx, ly, ran, ux, uy
INTEGER    :: i, j

!     EXTERNAL FUNCTIONS
! REAL :: ran
! EXTERNAL ran

!     COMMON BLOCKS
! COMMON /polyg/ nvs, vert, edges, np, totnvs

!     DEFINE BOX CONSTRAINTS FOR THE CENTRAL POINT

lx = 0.40D0 * (nx-1) * xstep
ux = 0.60D0 * (nx-1) * xstep
ly = 0.40D0 * (ny-1) * ystep
uy = 0.60D0 * (ny-1) * ystep

!     DEFINE CENTRAL POINT POLYGON

nvs(1) = 4

vert(1) = lx
vert(2) = ly
vert(3) = lx
vert(4) = uy
vert(5) = ux
vert(6) = uy
vert(7) = ux
vert(8) = ly

edges(1) = -1.0D0
edges(2) = 0.0D0
edges(3) = lx
edges(4) = 0.0D0
edges(5) = 1.0D0
edges(6) = -uy
edges(7) = 1.0D0
edges(8) = 0.0D0
edges(9) = -ux
edges(10) = 0.0D0
edges(11) = -1.0D0
edges(12) = ly

!     DEFINE CITY-POLYGONS CENTERED AT THE GRID POINTS AND OUTSIDE
!     THE CENTRAL REGION

np = 1
totnvs = 4

DO i = 0, nx - 1
  DO j = 0, ny - 1

    cx = i * xstep
    cy = j * ystep

    CALL RANDOM_NUMBER( ran )
    IF ((cx < lx-xstep .OR. cx > ux+xstep .OR. cy < ly-ystep .OR.  &
         cy > uy+ystep) .AND. ran <= prob) THEN

! GENERATE A NEW POLYGON

      np = np + 1
      CALL genpol(cx, cy, nvsvmi, nvsvma, rmin, rmax)
      totnvs = totnvs + nvs(np)

    END IF

  END DO
END DO

RETURN
END SUBROUTINE genpro



SUBROUTINE genpol(cx, cy, nvsvmi, nvsvma, rmin, rmax)

! N.B. Argument SEED has been removed.

!  This subroutine generates a polygon in R^2 with its vertices in a sphere
!  centered at point (cx,cy).  The number of vertices is randomly generated
!  satisfying nvsvmi <= number of vertices <= nvsvma.  The ratio of the
!  sphere is also randomly generated satisfying rmin <= ratio < rmax.
!  The generated polygon is stored in the common block "polyg".

!  On Entry:

!  cx    REAL (dp),
!        first coordinate of the center of the sphere,

!  cy    REAL (dp),
!        second coordinate of the center of the sphere,

!  nvsvmi integer,
!        minimum number of vertices,

!  nvsvma integer,
!        maximum number of vertices,

!  rmin double,
!        minimum ratio of the sphere,

!  rmax double,
!        maximum ratio of the sphere,

!  seed double,
!        seed for the random generation.

!  On Output:

!  The generated polygon is stored in the polyg common block described below.

!  Common block polyg:

!  common /polyg/nvs,vert,edges,np,totnvs

!  This structure represents, at any time, np polygons.
!  Position i of array nvs indicates the number of vertices of polygon i.
!  Arrays vert and edges store the vertices and edges of the polygons.

!  For example, if nvs(1) = 3 it indicates that the first polygon has
!  3 vertices (edges).  Then, if the vertices are (x1,y1), (x2,y2) and
!  (x3,y3), we have that vert(1) = x1, vert(2) = y1, vert(3) = x2,
!  vert(4) = y2, vert(5) = x3, and vert(6) = y3.  And, if the edges
!  (written as ax + by + c = 0) are a1 x + b1 y + c1 = 0,
!  a2 x + b2 y + c2 = 0, and a3 x + b3 y + c3 = 0 then edges(1) = a1,
!  edges(2) = b1, edges(3) = c1, edges(4) = a2, edges(5) = b2, edges(6) = c2,
!  edges(7) = a3, edges(8) = b3 and edges(9) = c3.

!  totnvs indicates the total number of vertices of the set of polygons.
!  This information is used when a new polygon is created to know the first
!  free position of arrays vert and edges at which the vertices and edges
!  of the new polygon will be saved.

!  Two additional details:

!  1) For each polygon, the vertices are ordered clockwise and edge i
!  corresponds to the edge between vertices i and i+1 (0 if i=n).

!  2) For each edge of the form ax + bx + c = 0, constants a, b and c are
!  chosen in such a way that (|a| = 1 or |b| = 1) and (a cx + b cy + c <= 0).

!     PARAMETERS

REAL (dp), INTENT(IN)  :: cx
REAL (dp), INTENT(IN)  :: cy
INTEGER, INTENT(IN)    :: nvsvmi
INTEGER, INTENT(IN)    :: nvsvma
REAL (dp), INTENT(IN)  :: rmin
REAL (dp), INTENT(IN)  :: rmax

REAL (dp), PARAMETER :: pi = 3.141592653589793_dp

!     COMMON SCALARS
! INTEGER :: np, totnvs

!     COMMON  ARRAYS
! REAL (dp)  :: edges(npmax*nvsmax*3), vert(npmax*nvsmax*2)
! INTEGER    :: nvs(npmax)

!     LOCAL SCALARS
REAL (dp)  :: r, ran
INTEGER    :: i

!     LOCAL ARRAYS
REAL (dp)  :: angl(nvsmax)

!     EXTERNAL FUNCTIONS
! REAL :: ran
! EXTERNAL ran

!     COMMON BLOCKS
! COMMON /polyg/ nvs, vert, edges, np, totnvs

!     GENERATE THE NUMBER OF VERTICES

CALL RANDOM_NUMBER( ran )
nvs(np) = nvsvmi + (nvsvma-nvsvmi+1)*ran

!     GENERATE THE RATIO OF THE SPHERE

CALL RANDOM_NUMBER( ran )
r = rmin + (rmax-rmin) * ran

!     GENERATE ALL ANGLES SATISFYING 0 <= ANGLE_I < 2*PI

DO i = 1, nvs(np)
  CALL RANDOM_NUMBER( ran )
  angl(i) = 2 * pi * ran
END DO

!     CLASSIFY THE ANGLES IN DECREASING ORDER

CALL class(nvs(np), angl)

!     CONSTRUCT THE VERTICES

DO i = 1, nvs(np)
  vert(2*(totnvs+i)-1) = cx + r * COS(angl(i))
  vert(2*(totnvs+i)) = cy + r * SIN(angl(i))
END DO

!     CONSTRUCT THE EDGES

DO i = totnvs + 1, totnvs + nvs(np) - 2
  CALL constr(vert(2*i-1), vert(2*i), vert(2*i+1), vert(2*i+2), vert(2*i+3), &
              vert(2*i+4), edges(3*i-2), edges(3*i-1), edges(3*i))
END DO

i = totnvs + nvs(np) - 1

CALL constr(vert(2*i-1), vert(2*i), vert(2*i+1), vert(2*i+2),  &
            vert(2*totnvs+1), vert(2*totnvs+2), edges(3*i-2), edges(3*i-1), &
            edges(3*i))

i = totnvs + nvs(np)

CALL constr(vert(2*i-1), vert(2*i), vert(2*totnvs+1), vert(2*totnvs+2),  &
            vert(2*totnvs+3), vert(2*totnvs+4), edges(3*i-2), edges(3*i-1), &
            edges(3*i))
RETURN
END SUBROUTINE genpol



SUBROUTINE constr(x1, y1, x2, y2, x3, y3, a, b, c)

!  This subroutine computes the real constants a, b and c of
!  the straight line ax + by + c = 0 in R^2 defined by the
!  points (x1,y1) and (x2,y2); such that the point (x3,y3)
!  satisfies the constraint ax + by + c <= 0.

!  On Entry:

!  x1    REAL (dp),
!        first coordinate of point (x1,y1),

!  y1    REAL (dp),
!        second coordinate of point (x1,y1),

!  x2    REAL (dp),
!        first coordinate of point (x2,y2),

!  y2    REAL (dp),
!        second coordinate of point (x2,y2),

!  x3    REAL (dp),
!        first coordinate of point (x3,y3),

!  y3    REAL (dp),
!        second coordinate of point (x3,y3).

!  On Return

!  a,b,c REAL (dp)
!        the desired constants.

REAL (dp), INTENT(IN)   :: x1
REAL (dp), INTENT(IN)   :: y1
REAL (dp), INTENT(IN)   :: x2
REAL (dp), INTENT(IN)   :: y2
REAL (dp), INTENT(IN)   :: x3
REAL (dp), INTENT(IN)   :: y3
REAL (dp), INTENT(OUT)  :: a
REAL (dp), INTENT(OUT)  :: b
REAL (dp), INTENT(OUT)  :: c

IF (x1 == x2 .AND. y1 == y2) THEN
  WRITE (*,FMT=*) 'ERROR IN FUNCTION CONSTRAINT: X1=X2 AND Y1=Y2'
END IF

IF (y1 /= y2) THEN
  a = 1.0D0
  b = -(x2-x1) / (y2-y1)
  c = -(x1 + b*y1)

ELSE
  a = 0.0D0
  b = 1.0D0
  c = -y1
END IF

IF (a*x3 + b*y3 + c > 0.0D0) THEN
  a = -a
  b = -b
  c = -c
END IF

RETURN
END SUBROUTINE constr



SUBROUTINE class(n,x)

!  This subroutine classifies the elements of a vector in decreasing order,
!  i.e., on output: x(1) >= x(2) >= ... >= x(n).

!  On Entry:

!  n     integer,
!        number of elements of the vector to be classified,

!  x     REAL (dp) x(n),
!        vector to be classified.

!  On Return

!  x     REAL (dp) x(n),
!        classified vector.

!     ARGUMENTS

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(n)

!     LOCAL SCALARS
REAL (dp)  :: aux, xmax
INTEGER    :: i, j, pos

DO i = 1, n

  xmax = x(i)
  pos = i

  DO j = i + 1, n
    IF (x(j) > xmax) THEN
      xmax = x(j)
      pos = j
    END IF
  END DO

  IF (pos /= i) THEN
    aux = x(i)
    x(i) = x(pos)
    x(pos) = aux
  END IF

END DO

RETURN
END SUBROUTINE class

END MODULE Common_SPGMA2



PROGRAM spgma2
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-05-04  Time: 17:43:54

!  SPG Test driver. This master file uses SPG to solve the location problem
!  described in the Test Problems section.
!  It generates the problems, calls the optimizer to solve them and writes
!  the reports (tables).

!  This version 17 JAN 2000 by E.G.Birgin, J.M.Martinez and M.Raydan.
!  Reformatted 03 OCT 2000 by Tim Hopkins.
!  Final revision 03 JUL 2001 by E.G.Birgin, J.M.Martinez and M.Raydan.

USE Spectral_Projected_Grad
USE Common_SPGMA2
IMPLICIT NONE

!     LOCAL SCALARS
REAL (dp)  :: pginfn, pgtwon, eps, eps2, f, prob, rmax, rmin, xstep, ystep
REAL       :: time
INTEGER    :: fcnt, flag, gcnt, iter, m, maxfc, maxit, n, nvsvma,  &
              nvsvmi, nx, ny, pnum
LOGICAL    :: output

!     LOCAL ARRAYS
REAL (dp)  :: x(nmax)
REAL       :: dum(2)

!     EXTERNAL SUBROUTINES
! EXTERNAL dtime, genpro, spg

!     DEFINE THE PROBLEM

OPEN(UNIT=8, FILE='DATA2', STATUS='OLD')

DO

   WRITE (*, FMT=5000)
   READ (8, FMT=*, END=20) pnum, nx, ny, prob, nvsvmi, nvsvma

   IF (nvsvmi < 3) THEN
     WRITE (*, FMT=*) 'NVSVMI MUST BE GREATER THAN OR EQUAL TO 3'
     READ (*, FMT=*)
   END IF

   IF (nvsvma > nvsmax) THEN
     WRITE (*, FMT=*) 'NVSVMA MUST BE LESS THAN OR EQUAL TO', nvsmax
     WRITE (*, FMT=*) 'REDUCE NVSVMA OR INCREASE NVSMAX'
     READ (*, FMT=*)
   END IF

   IF (nx*ny*prob > npmax) THEN
     WRITE (*, FMT=*) 'NX*NY*PROB MUST BE LESS THAN ', npmax
     WRITE (*, FMT=*) 'REDUCE NX, NY OR PROB OR INCREASE NPMAX'
     READ (*, FMT=*)
   END IF

   xstep = 5.0D0
   ystep = 5.0D0

   rmin = 1.0D0
   rmax = 2.0D0

   !     GENERATE THE PROBLEM

   CALL genpro(nx, ny, xstep, ystep, prob, nvsvmi, nvsvma, rmin, rmax)

   WRITE (*, FMT=*) 'NUMBER OF POLYGONS: ', np
   WRITE (*, FMT=*) 'NUMBER OF VERTICES (AND EDGES): ', totnvs

   !     DEFINE INITIAL POINT

   n = 2 * np
   x(1:n) = 0.0_dp

   !     SET UP THE INPUT DATA OF THE OPTIMIZATION ALGORITHM

   output = .false.
   maxit = 1000
   maxfc = 2000
   eps = 0.0D0
   eps2 = 1.0D-06
   m = 10

   !     CALL THE OPTIMIZER

   CALL CPU_TIME(dum(1))
   CALL spg(n, x, m, eps, eps2, maxit, maxfc, output, f, pginfn, pgtwon,  &
            iter, fcnt, gcnt, flag)
   CALL CPU_TIME(dum(2))

   time = dum(2) - dum(1)

   !     WRITE STATISTICS

   WRITE (*,FMT=5100) f, pginfn, SQRT(pgtwon), flag
   WRITE (*,FMT=5200) iter, fcnt, gcnt
   WRITE (*,FMT=5300) time

   !     WRITE SOLUTION

END DO

20 STOP

5000 FORMAT (/'---------------------------------------'//)
5100 FORMAT (/' F = ', g17.10 / ' PGINFNORM = ', g16.10 /  &
             ' PGTWONORM = ', g16.10 / ' FLAG = ', i1)
5200 FORMAT (/' ITER = ', i10,/' FCNT = ', i10 / ' GCNT = ', i10)
5300 FORMAT (/' TIME = ', f12.2, ' Secs.')
END PROGRAM spgma2



SUBROUTINE evalf(n, x, f, inform)

!  This subroutine computes the objective function.

!  On Entry:

!  n     integer,
!        size of the problem,

!  x     REAL (dp) x(n),
!        point at which the function will be evaluated.

!  On Return

!  f     REAL (dp),
!        function value at x,

!  inform integer,
!        termination parameter:
!        0 = the function was successfully evaluated,
!        1 = some error occurs in the function evaluation.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(n)
REAL (dp), INTENT(OUT)  :: f
INTEGER, INTENT(OUT)    :: inform

!     LOCAL SCALARS
REAL (dp)  :: diff1, diff2, dist
INTEGER    :: i, ndist

inform = 0

f = 0.0_dp
ndist = n / 2 - 1
DO i = 1, ndist
  diff1 = x(1) - x(2*i+1)
  diff2 = x(2) - x(2*i+2)
  dist = SQRT(diff1**2 + diff2**2)
  IF (dist <= 1.0D-4) THEN
    WRITE (*,FMT=*) 'ERROR IN PROBLEM DEFINITION (DIST TOO SMALL)'
    inform = 1
  END IF
  f = f + dist
END DO

f = f / ndist
RETURN
END SUBROUTINE evalf



SUBROUTINE evalg(n, x, g, inform)

!  This subroutine computes the gradient of the objective function.

!  On Entry:

!  n     integer,
!        size of the problem,

!  x     REAL (dp) x(n),
!        point at which the gradient will be evaluated.

!  On Return

!  g     REAL (dp) g(n),
!        gradient vector at x,

!  inform integer,
!        termination parameter:
!        0 = the gradient was successfully evaluated,
!        1 = some error occurs in the gradient evaluation.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(n)
REAL (dp), INTENT(OUT)  :: g(n)
INTEGER, INTENT(OUT)    :: inform

!     LOCAL SCALARS
REAL (dp)  :: diff1, diff2, dist
INTEGER    :: i, ndist

inform = 0

g(1) = 0.0D0
g(2) = 0.0D0
ndist = n / 2 - 1
DO i = 1, ndist
  diff1 = x(1) - x(2*i+1)
  diff2 = x(2) - x(2*i+2)
  dist = SQRT(diff1**2 + diff2**2)
  g(2*i+1) = -(diff1/dist) / ndist
  g(2*i+2) = -(diff2/dist) / ndist
  g(1) = g(1) - g(2*i+1)
  g(2) = g(2) - g(2*i+2)
END DO
RETURN
END SUBROUTINE evalg



SUBROUTINE proj(n, x, inform)

!  This subroutine computes the projection of an arbitrary point onto the
!  feasible set.  Since the feasible set can be described in many ways,
!  its information is in a common block.

!  In this particular implementation of the projection subroutine for the
!  location problem, each point z_i in R^2 (stored at positions 2i-1 and 2i
!  of x) must be projected onto polygon P_i.
!  See the Test Problems section for details.

!  On Entry:

!  n     integer,
!        size of the problem,

!  x     REAL (dp) x(n),
!        point that will be projected.

!  On Return

!  x     REAL (dp) x(n),
!        projected point,

!  inform integer,
!        termination parameter:
!        0 = the projection was successfully done,
!        1 = some error occurs in the projection.

!     PARAMETERS

USE Common_SPGMA2
IMPLICIT NONE

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(n)
INTEGER, INTENT(OUT)       :: inform

!     COMMON ARRAYS
! REAL (dp)  :: edges(npmax*nvsmax*3), vert(npmax*nvsmax*2)
! INTEGER    :: nvs(npmax)

!     LOCAL SCALARS
REAL (dp)  :: xproj, yproj
INTEGER    :: base, i

!     COMMON BLOCKS
! COMMON /polyg/ nvs, vert, edges, np, totnvs

inform = 0

base = 0

DO i = 1, np
  
! PROJECT z_i ONTO P_i
  
  CALL projp(i, base, x(2*i-1), x(2*i), xproj, yproj, inform)
  
  x(2*i-1) = xproj
  x(2*i) = yproj
  
  base = base + nvs(i)
  
END DO
RETURN
END SUBROUTINE proj
