c**********************************fill_land.f********************************
c This program takes our product of Levitus formatted, land-masked data and
c fills in the land values using a poisson solver (in the extrap.f subroutine
c taken directly from the MOM2 source code).  The final data are written to
c a Levitus formatted file of the same name with a .fil appendage, which can 
c then be used with a model grid of any resolution.
c To use this code, you must:
c
c 1. Find this line near the top of the code:
c	 parameter (idim=360,jdim=180,kdim=33,ct=1.e-5,
c     *              bad=-99.,value=2,ms=500)
c 2. In this line, set "value" equal to 1 for temp or 2 for salt. 
c 3. In this line, set "kdim" to either 24 (monthly) or 33 (seasonal and annual)
c 3. Find this line in the code: data month/ '00'/
c    Change the '00' to reflect the appropriate season or month. 
c    ('00' is annual data, '13' is winter, '15' is summer, '01' through '12'
c     are months)
c 4. Make sure the filename is appropriate for your data. It will be if you have
c    not changed the name of the files, and this program is in the same directory
c    with them. If you change filename, make sure the filename
c     variable is long enough to accomodate the full path name: 
c     e.g:character*90 filename.
c
c    Currently, filename=title//month//'_p3.obj' 
c       where title='Salt' or 'Temp' and month is set by you as explained in
c       #3 above.
c 4. Compile and run the program.
c    >f77 fill_land.f (to compile)
c    >a.out           (to run)

c***************************************************************************
c temp: character*15 filename,character*14 file1,character*11 dir,character*4 title;
c salt:character*8 filename,character*12 file1,character*8 dir,character*3 title


c ====EDIT parameter TO SPECIFY WHAT TYPE OF DATA AND # OF LEVELS====
c --> value =1 for temp or 2 for salt
c --> kdim = 33 for annual, winter or summer fields
c --> kdim = 24 for monthly fields

	 parameter (idim=360,jdim=180,kdim=24,ct=1.e-5,
     *              bad=-99.,value=2,ms=500)

         character*15 filename
         character*19 file1
         character*11 dir
         character*2 month
         character*4 suff
         character*4 title
       	 real tmp(idim,jdim),sor(idim,jdim),
     *        res(idim,jdim),jnk(idim,jdim)
	 integer gt,init,ilevmsk(idim,jdim)
	 
	 
	 gt=1      	

c =======EDIT THIS LINE TO CHOOSE THE DATA FILE YOU WISH TO FILL======
c  ('00' is annual data. '01' is Jan. '13' is winter. '15' is summer)
         
        data month/ '01'/
c=====================================================================    


c temperature or salinity?

	if (value .eq. 1) then
	 dir='Temperature'
	 title='Temp'
	 suff='temp'
	else
	 dir='Salinity'
	 title='Salt'
	 suff='salt'
	endif
            	 

c for the very first depth, set tmp over land to be equal to
c zero.  For every other depths, use the previous depth's value
c to initialize extrap.
	init=1
	

c Input OI data file. Specify the path name to your directory
c ====YOU MAY NEED TO CHANGE THE FILENAME. THE LINE BELOW SHOULD WORK
c ====IF YOU HAVE NOT CHANGED THE NAMES OF THE FILES. FILES NEED TO BE
c ====IN THE SAME DIRECTORY AS THE PROGRAM. THE OUTPUT WILL NOT OVERWRITE
C ====THE ORIGINAL DATA. A NEW FILE WITH THE SAME NAME WILL BE CREATED
C ====WITH A .fil APPENDAGE.

 	 filename=title//month//'_p3.obj'
     
	write(*,*) filename
         open(22,file=filename,status='old',form='formatted') 
              

 	 file1=filename//'.fil'
	  write(*,*) file1
	  open(23,file=file1,
     *		form='formatted')
     
        do 50 k=1,kdim
       

c initialize landmask to water for each depth 
	do 24 i=1,idim
	 do 23 j=1,jdim
	    ilevmsk(i,j)=1
23	continue
24	continue	    

           print *, jnk(i,j), jdim 
           read(22,80) ((jnk(i,j),i=1,360),j=1,jdim)        
   

c missing point: if surface, fill tmp with 0, otherwise fill with previous
c non-missing depth value. ilevmsk is appropriately filled for any depth.
	   do 31 i=1,idim
	   do 30 j=1,jdim
	    if (jnk(i,j) .le. bad) then
	    ilevmsk(i,j)=0
	     if (init .eq. 1) then
	      tmp(i,j)=0.0
	     endif
	    else
	    tmp(i,j)=jnk(i,j)
	    endif
30	   continue
31	   continue

	  call extrap(tmp,ilevmsk,sor,res,idim,jdim,ms,ct,
     *		      'levitus data',gt)

	  

     	
	write(23,80) tmp
	init=init+1
50	 continue	   	  

80     format(10f8.4)       	
       	   
	close(22)	
	close(23)

	 end
	 

c This routine is taken directly from 2MOM2.f and is used there as well to
c fill in land values using a poisson solver, so that the data can be
c processed by the MOM2 model code.
c*************************************************************************
c The extrap.f subroutine is taken from the MOM source code in
c util.f.  It has been modified in this way:
c The stdunits.h library is not used, and stdout is replaced with file
c    unit 6
c The if statment for cyclic conditions has also been removed.
c
c

      subroutine extrap (a, land, sor, res, il, jl, maxscn, crit, text
     &,                  gtype)
c
c=======================================================================
c
c     utility to extrapolate values into land areas neglecting
c     non-uniformity or asymmetry in the grid by solving a simple
c     heat eqn: del**2(a) = 0 over land areas using values over ocean
c     areas as boundary conditions.
c     this alleviates the problem of mismatched land/sea areas due to
c     different geometries or resolutions when interpolating between
c     atmospheric and ocean model grids.
c     the intent is to force reasonable values into land areas
c     near coastlines. far from coasts, the extrapolations may not be
c     reasonable.
c
c     note: the values over land are used as an initial guess field
c           and need to be specified
c
c
c     inputs:
c
c     a       = array with land areas to be filled. land areas contain
c               initial guess field.
c     land    = mask = (0, non zero) to indicate (land, non land) area
c     il      = number of points along 1st dimension to be filled
c     jl      = number of points along 2nd dimension to be filled
c     maxscn  = maximum number of passes allowed in relaxation
c     crit    = criterion for ending relaxation before "maxscn" limit
c     text    = character string (up to 15 chars) to identify data
c     gtype   = grid type = (1,2) to identify (ocean, atmosphere) grid
c     sor     = scratch area
c     res     = scratch area
c
c
c     outputs:
c
c     a       = array with extrapolated values in land areas.
c               non land areas remain unchanged.
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      logical done
c#include "stdunits.h"
      integer gtype
      character*(*) text
      parameter (c0=0.0, p25=0.25)
      dimension a(il,jl), land(il,jl), res(il,jl), sor(il,jl)



c
c-----------------------------------------------------------------------
c
c     solve a simple poisson eqn by relaxation to extrapolate data into
c     land areas using values over non land areas as boundary values.
c
c     note: sucessive calls to extrap will require fewer scans beacuse
c           the initial guess field over land areas gets better with
c           each call.
c-----------------------------------------------------------------------
c
c
c     check on the grid type: atmosphere or ocean
c
      if (gtype .ne. 1 .and. gtype .ne. 2) then
        write (6,98) gtype
        stop '=>extrap'
      endif     
c
c-----------------------------------------------------------------------
c     set the relaxation coefficient to zero over ocean or air
c     relc is somewhat arbitrary
c-----------------------------------------------------------------------
c
      relc = 0.6
      do j=1,jl
        do i=1,il
          if (land(i,j) .eq. 0) then
            sor(i,j) = relc
          else
            sor(i,j) = c0
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     iterate until errors are acceptable.
c-----------------------------------------------------------------------
c     
      n = 0
100   continue
        resmax = c0
        done   = .true.
        n    = n + 1
        do j=2,jl-1
          do i=2,il-1
            res(i,j) = p25*(a(i-1,j) + a(i+1,j) + a(i,j-1) + a(i,j+1)) 
     &                 - a(i,j)
          enddo
        enddo
        do j=2,jl-1
          do i=2,il-1
            res(i,j) = res(i,j)*sor(i,j)
            a(i,j) = a(i,j) + res(i,j)
            absres = abs(res(i,j))
            if (absres .gt. crit) done = .false.
            resmax = max(absres,resmax)
          enddo
        enddo
c
c-----------------------------------------------------------------------
c       set conditions at edge of grid
c-----------------------------------------------------------------------
c
        if (gtype .eq. 1) then
c
c         use cyclic or no flux conditions on ocean grids
c
          do j=1,jl
            a(1,j)  = a(il-1,j)
            a(il,j) = a(2,j)

          enddo
        elseif (gtype .eq. 2) then
c
c         always put cyclic conditions on atmosphere grids
c
          do j=1,jl
            a(1,j)  = a(il-1,j)
            a(il,j) = a(2,j)
          enddo
        endif
c
c       no flux condition at northern and southern boundaries
c
        do i=1,il
          a(i,1)  = a(i,2)
          a(i,jl) = a(i,jl-1)
          enddo
c
      if (.not. done .and. n .le. maxscn) go to 100
c
      write (6,99) text, n, resmax

   
99    format (1x,'==> Extrapolated ',a15,' into land using ',i4
     &,       ' scans.  max residual=', g14.7)
98    format (1x,'==> Error:   gtype =',i6,' in extrap')

97	format (g14.7)

      return
      
      end


	  
	 
