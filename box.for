*                       File : box.for

*------------------- Box -----------------------------

        program MAIN
         include 'box.inc'
         integer fw, num, ng, set_read, p
         character*12 FileSpecIn, FileSpecOut
         character*12 FileFldIn, FileFldOut

          write(*,*) ' Very Slow Fourier Transform (VSFT)'
          fw = 12
          call OpenInputFile( fw, 'box.in','formatted')
          read( fw, * ) set_read
          read( fw, * ) ng
          read( fw, * ) p
          read( fw, * )
          read( fw, * ) NumK
          read( fw,'(a12)') FileSpecIn
          read( fw,'(a12)') FileFldOut
          read( fw, * )
          read( fw, * ) NumKout
          read( fw,'(a12)') FileSpecOut
          read( fw,'(a12)') FileFldIn
          close(fw)
          
          write(*,*) ng*(p+1), ng*p+1

          if(set_read.Eq.0) then
            ng = ng*(p+1)
          else
            ng = ng*p+1
          endif

          if((NumK.Gt.MaxK).Or.(NumKout.Gt.MaxK)) then
            write(*,'(''Error, max K ='',i5)') MaxK
            stop
          endif

C           Pi_const = 3.14159265358979323846D0
          Pi_const = 3.141592653589793115998D0
          call GenerGrid( ng, p , set_read)

          if( set_read .Eq. 1 ) then
            call ReadFld( FileFldIn )
            write(*,*) ' Extract spectrum, please wait. '
            call Spectrum( )
            call WriteSpectrum( FileSpecOut )
          else
            call ReadSpectrum( FileSpecIn, num )
            write(*,*) ' Generation velocity fields, please wait '
            call InitFourier( num )
            call GenerVel( )
            call WriteFld( ng, p, FileFldOut)
          endif
        end



*-----------------------  Generation ------------------------------

        subroutine VelFun( mk, vk, ik, lf, iran )
         include 'box.inc'
         real *16 mk, vk(*)
         integer ik, lf, iran
         integer i, j, k
         real *16 dv(3), arg, ang, fcos
         real *16 Random, norm
         integer n

          ang = Random(iran)*2.0*Pi_const
          do n = 1, Ngrd, 1
            dv(n) = Div(lf,n)
          enddo
          norm = sqrt(SpecE(ik)/(2.0*CtK(ik)))

          do k = Nk1, Nk2, 1
          do j = Nj1, Nj2, 1
          do i = Ni1, Ni2, 1
            arg  = vk(1)*Xg(i,j,k) + vk(2)*Yg(i,j,k) + vk(3)*Zg(i,j,k)
            fcos = 2.0 * norm * cos(arg+ang)
            do n = 1, Ngrd, 1
              Vel(i,j,k,n) = Vel(i,j,k,n) + fcos * dv(n)
            enddo
          enddo
          enddo
          enddo
        end

        subroutine AmpFun( mk, vk, ik, lf, iran )
         include 'box.inc'
         real *16 mk, vk(*)
         integer ik, lf, iran
         integer n, type
         real *16 Random, GaussRandom

          type = 1
          do n = 1, Ngrd, 1
            if( type .Eq. 1 ) then
              Amp(lf,n) = GaussRandom(iran)
            else if( type .Eq. 2 ) then
              Amp(lf,n) = 2.0*Random(iran)
            else
              Amp(lf,n) = 1.0
            endif
          enddo
        end

        subroutine DivFun( mk, vk, ik, lf, iran )
         include 'box.inc'
         real *16 mk, vk(*)
         integer ik, lf, iran
         real *16 dd, sum, am(3)
         integer n, m

          do n = 1, Ngrd, 1
            am(n) = Amp(lf,n)
          enddo
          do n = 1, Ngrd, 1
            sum = 0.0
            do m = 1, Ngrd, 1
              dd = 0.0
              if(m .Eq. n) dd = 1.0
              sum = sum + (dd - vk(n)*vk(m)/mk**2)*am(m)
            enddo
            Div(lf,n) = sum
          enddo
        end

        subroutine AmpNormFun( mk, vk, ik, lf, iran )
         include 'box.inc'
         real *16 mk, vk(*)
         integer ik, lf, iran
         real *16 tmp
         integer n

          do n = 1, Ngrd, 1
            tmp = sqrt(1.0D0/AmpNorm(ik,n))
            Amp(lf,n) = Amp(lf,n) * tmp
          enddo
        end

        subroutine DivNormFun( mk, vk, ik, lf, iran )
         include 'box.inc'
         real *16 mk, vk(*)
         integer ik, lf, iran
         real *16 tmp
         integer n

          tmp = sqrt(2.0D0/DivNorm(ik))
          do n = 1, Ngrd, 1
            Div(lf,n) = Div(lf,n) * tmp
          enddo
        end


        subroutine CountFun( mk, vk, ik, lf, iran )
         include 'box.inc'
         real *16 mk, vk(*)
         integer ik, lf, iran
         real *16 sum
         integer n

          sum = 0.0
          do n = 1, Ngrd, 1
            sum = sum + Div(lf,n)**2
          enddo
          if(sum .Gt. 1.0D-10) then
            do n = 1, Ngrd, 1
              AmpNorm(ik,n) = AmpNorm(ik,n) + Amp(lf,n)**2
            enddo
            DivNorm(ik) = DivNorm(ik) + sum
            CtK(ik) = CtK(ik) + 1.0
          endif
        end

        subroutine Count(  )
         include 'box.inc'
         integer k, n
         external CountFun

          do k = 0, MaxK, 1
            do n = 1, Ngrd, 1
              AmpNorm(k,n) = 0.0
            enddo
            DivNorm(k) = 0.0
            CtK(k) = 0.0
          enddo
          call FourierCycle( CountFun, NumK )
          do k = 1, NumK, 1
            DivNorm(k) = DivNorm(k)/CtK(k)
            do n = 1, Ngrd, 1
              AmpNorm(k,n) = AmpNorm(k,n)/CtK(k)
            enddo
          enddo
        end

        subroutine GenerVel(  )
         include 'box.inc'
         integer i, j, k, n, lf
         integer size
         external DivFun, AmpFun
         external DivNormFun, AmpNormFun
         external VelFun

          size = (2*NumK+1)**3
          do lf = 0, size-1, 1
            do n = 1, Ngrd, 1
              Div(lf,n) = 0.0
              Amp(lf,n) = 0.0
            enddo
          enddo

          do k = Nk1, Nk2, 1
          do j = Nj1, Nj2, 1
          do i = Ni1, Ni2, 1
            do n = 1, Ngrd, 1
              Vel(i,j,k,n) = 0.0
            enddo
          enddo
          enddo
          enddo

          call FourierCycle( AmpFun, NumK )
          call FourierCycle( DivFun, NumK )
          call Count( )
          call FourierCycle( AmpNormFun, NumK )
          call FourierCycle( DivFun, NumK )
          call Count( )
          call FourierCycle( DivNormFun, NumK )

          call FourierCycle( VelFun, NumK )
        end

        subroutine InitFourier( num )
         include 'box.inc'
         integer num, k, n1, n2, n
         real *16 rk, rk1, rk2, e1, e2, sum

          if(Kin(1).Gt.0.5D0) then
            write(*,*)  'Error: incorrect energy profiles'
            stop
          endif
          SpecE(0) = 0.0
          do k = 1, NumK, 1
            rk  = k
            rk1 = rk - 0.5D0
            rk2 = rk + 0.5D0
            n1 = 1
            n2 = 1
            do while(Kin(n1).Lt.rk1)
              n1 = n1 + 1
            enddo
            do while(Kin(n2).Lt.rk2)
              n2 = n2 + 1
            enddo

            e1 = Ein(n1) + (Ein(n1+1)-Ein(n1))*
     &           (rk1-Kin(n1))/(Kin(n1+1)-Kin(n1))

            e2 = Ein(n2-1) + (Ein(n2)-Ein(n2-1))*
     &           (rk2-Kin(n2-1))/(Kin(n2)-Kin(n2-1))

            sum = 0.5*(e1+Ein(n1))*(Kin(n1)-rk1) +
     &            0.5*(e2+Ein(n2-1))*(rk2-Kin(n2-1))

            do n = n1+1, n2-1, 1
              sum = sum + 0.5*(Ein(n)+Ein(n-1))*(Kin(n)-Kin(n-1))
            enddo
            SpecE(k) = sum
          enddo
        end


*-----------------------  Spectrum ------------------------------

        subroutine SpectrumFun( mk, vk, ik, lf, iran )
         include 'box.inc'
         real *16 mk, vk(*)
         integer i, j, k
         integer ik, lf, iran, n
         real *16 sum_r(3), sum_i(3)
         real *16 vl, tmp, fcos, fsin, arg

          do n = 1, Ngrd, 1
            sum_r(n) = 0.0
            sum_i(n) = 0.0
          enddo
          do k = Nk1+1, Nk2, 1
          do j = Nj1+1, Nj2, 1
          do i = Ni1+1, Ni2, 1
            arg = vk(1)*Xg(i,j,k)+vk(2)*Yg(i,j,k)+vk(3)*Zg(i,j,k)
            fcos = cos(-arg)
            fsin = sin(-arg)
            do n = 1, Ngrd, 1
              vl = Vel(i,j,k,n)
              sum_r(n) = sum_r(n) + vl * fcos
              sum_i(n) = sum_i(n) + vl * fsin
            enddo
          enddo
          enddo
          enddo
          tmp = 0.0D0
          do n = 1, Ngrd, 1
            tmp = tmp + sum_r(n)**2 + sum_i(n)**2
          enddo
          EnK(ik) = EnK(ik) + 2.0*tmp
        end

        subroutine Spectrum( )
         include 'box.inc'
         integer k
         real *16 tmp
         external SpectrumFun

          do k = 0, NumKout, 1
            EnK(k) = 0.0
          enddo
          call FourierCycle(SpectrumFun, NumKout)

          tmp = 0.5D0/VolBox**2
          do k = 1, NumKout, 1
            EnK(k) = EnK(k) * tmp
          enddo
        end

*--------------------  Utility  ----------------------------------

        subroutine FourierCycle( FourierFun, num_f )
         integer num_f
         integer nf1, nf2, nf3
         integer szf, lf, iran, ik
         real *16 mk_max, mk, vk(3), rk, eps
         external FourierFun

          eps = 1.0D-15
          mk_max = dble(num_f) + 0.5D0
          szf    = 2*num_f + 1
          iran = 1
          do nf1 = 0, num_f, 1
c            write(*, *) nf1
            vk(1) = nf1
          do nf2 = -num_f, num_f, 1
            vk(2) = nf2
          do nf3 = -num_f, num_f, 1
            vk(3) = nf3
            if((nf1.Eq.0).And.( (nf2.Lt.0).Or.
     &          ((nf2.Eq.0).And.(nf3.Lt.0)) )) then
            else
              mk = sqrt(vk(1)**2+vk(2)**2+vk(3)**2)
              if((mk.Le.mk_max).And.(mk.Gt.eps)) then
                lf = nf1+num_f+(nf2+num_f+(nf3+num_f)*szf)*szf
                ik = 0
                rk = ik
                do while(.Not.((mk.Gt.rk-0.5D0).And.(mk.Le.rk+0.5D0)))
                 ik = ik + 1
                 rk = ik
                enddo
                call FourierFun( mk, vk, ik, lf, iran )
              endif
            endif
          enddo
          enddo
          enddo
        end


        function Random( code )
         real *16 Random
         integer code
         integer n1, n2, n3
         real *16 rm

          n1 = 714025
          n2 = 1366
          n3 = 150899
          rm = 1.0D0/n1
          code   = mod(n2*code+n3,n1)
          Random = code * rm
          code   = mod(n2*code+n3,n1)
        end

        function GaussRandom( code )
         real *16 GaussRandom, Random
         integer code
         real *16 rr, v1, v2

          rr = 2.0
          do while( rr .Gt. 1.0 )
            v1 = 2.0*Random(code)-1.0
            v2 = 2.0*Random(code)-1.0
            rr = v1**2 + v2**2
          enddo
          GaussRandom = v1 * sqrt(-2.0*log(rr)/rr)
        end

        subroutine OpenInputFile( fw, FileName, str )
         integer fw, ioStatus
         character*(*) FileName, str

          open( fw, file=FileName, status='OLD',
     ,          form = str, iostat=ioStatus )
          if( ioStatus .Ne. 0 ) then
            write( *, '(''# Cant open file : '',A12)') FileName
            stop
          else
            write( *, '(''# Open file : '',A12)') FileName
          endif
        end

        subroutine ReadFld( FileName )
         include 'box.inc'
         integer fw, i, j, k, n
         character*(*) FileName

          fw = 12
          call OpenInputFile( fw, FileName,'formatted')
          do k = Nk1, Nk2, 1
          do j = Nj1, Nj2, 1
          do i = Ni1, Ni2, 1
            read(fw,"(F19.16)") (Vel(i,j,k,n), n=1,Ngrd)
            !write(*,*) Vel(i,j,k,1),Vel(i,j,k,2),Vel(i,j,k,3)
          enddo
          enddo
          enddo
          close( fw )
        end

        subroutine WriteFld(ng, p, FileName )
         include 'box.inc'
         integer fw
         integer ng, i, j, k, l, m, n, ne, pp, p, set_read
         integer ll, mm, nn, o
         integer Ni3, Nj3, Nk3, indi, indj, indk
         real *16 dh, ri rj, rk, rl, rm, rn, rp
         real *16 X1, Y1, Z1, X2, Y2, Z2, wx, wy, wz
         real *16 Xo, Yo, Zo, U2, V2, W2
         character*(*) FileName

         !CONVERT FROM CARTESIAN GRID TO THE SOLUTION GRID
         ne = ng/(p+1)
         pp = p+1
         Ni1 = 1
         Nj1 = 1
         Nk1 = 1
         Ni2 = ng
         Nj2 = ng
         Nk2 = ng
         Ni3 = ne
         Nj3 = ne
         Nk3 = ne

         dh = (2.0 * Pi_const)/dble(ne)

         do k = Nk1, Nk3, 1
         do j = Nj1, Nj3, 1
         do i = Ni1, Ni3, 1
           do l = 1, pp, 1
           do m = 1, pp, 1
           do n = 1, pp, 1

             indi = (i-1)*pp+l
             indj = (j-1)*pp+m
             indk = (k-1)*pp+n

             !RESET INTERP VECTOR
             Vel2(indi,indj,indk,1) = 0
             Vel2(indi,indj,indk,2) = 0
             Vel2(indi,indj,indk,3) = 0

             !POSITION IN THE SOLUTION SPACE
             X1 = Xs(indi,indj,indk)
             Y1 = Ys(indi,indj,indk)
             Z1 = Zs(indi,indj,indk)

             !BUILD THE LAGRANGE INERPOLATION
             do ll = 1, pp, 1
             do mm = 1, pp, 1
             do nn = 1, pp, 1

             indi = (i-1)*pp+ll
             indj = (j-1)*pp+mm
             indk = (k-1)*pp+nn

             !POSITION IN THE CARTESIAN SPACE
             X2 = Xg(indi,indj,indk)
             Y2 = Yg(indi,indj,indk)
             Z2 = Zg(indi,indj,indk)

             U2 = Vel(indi,indj,indk,1)
             V2 = Vel(indi,indj,indk,2)
             W2 = Vel(indi,indj,indk,3)

             !BUILD THE WEIGHTS
             wx = 1
             wy = 1
             wz = 1
             
             do o = 1, pp, 1
               indi = (i-1)*pp+o
               indj = (j-1)*pp+mm
               indk = (k-1)*pp+nn
               
               Xo = Xg(indi,indj,indk)
               if(o.ne.ll) wx = wx*(X1-Xo)/(X2-Xo)
             enddo

             do o = 1, pp, 1
               indi = (i-1)*pp+ll
               indj = (j-1)*pp+o
               indk = (k-1)*pp+nn
               
               Yo = Yg(indi,indj,indk)
               if(o.ne.mm) wy = wy*(Y1-Yo)/(Y2-Yo)
             enddo

             do o = 1, pp, 1
               indi = (i-1)*pp+ll
               indj = (j-1)*pp+mm
               indk = (k-1)*pp+o
               
               Zo = Zg(indi,indj,indk)
               if(o.ne.nn) wz = wz*(Z1-Zo)/(Z2-Zo)
             enddo

             !ADD THEM UP
             indi = (i-1)*pp+l
             indj = (j-1)*pp+m
             indk = (k-1)*pp+n

             Vel2(indi,indj,indk,1) = Vel2(indi,indj,indk,1) +
     &            wx*wy*wz*U2
             Vel2(indi,indj,indk,2) = Vel2(indi,indj,indk,2) +
     &            wx*wy*wz*V2
             Vel2(indi,indj,indk,3) = Vel2(indi,indj,indk,3) +
     &            wx*wy*wz*W2
               
             enddo
             enddo
             enddo

           enddo
           enddo
           enddo
         enddo
         enddo
         enddo

          fw = 12
          open (unit=fw,file=Filename,action="write",status="replace")
          do k = Nk1, Nk2, 1
          do j = Nj1, Nj2, 1
          do i = Ni1, Ni2, 1
            write(fw,"(14F22.18)") Xs(i,j,k),Ys(i,j,k),
     &          Zs(i,j,k),(Vel2(i,j,k,n), n=1,Ngrd)
          enddo
          enddo
          enddo
          close( fw )
        end

        subroutine ReadSpectrum( FileName, num  )
         include 'box.inc'
         integer fw, i, num
         character*(*) FileName

          fw = 12
          call OpenInputFile( fw, FileName,'formatted')
          read(fw, * ) num
          if( num .Gt. MaxIn ) then
            write(*,*) ' Error: incorrect number energy profiles'
            stop
          endif
          do i = 1, num,  1
            read(fw,*) Kin(i), Ein(i)
          enddo
          close(fw)
        end

        subroutine WriteSpectrum( FileName  )
         include 'box.inc'
         integer fw, k
         real *16 rk
         character*(*) FileName

          fw = 12
          open( fw, file = FileName, form = 'formatted')
          do k = 1, NumKout,  1
            rk = k
            write(fw,'(1p2e14.5)') rk, EnK(k)
          enddo
          close(fw)
        end

        subroutine GenerGrid( ng, p, set_read)
         include 'box.inc'
         integer ng, i, j, k, l, m, n, ne, pp, p, set_read
         integer Ni3, Nj3, Nk3, indi, indj, indk
         real *16 dh, ri rj, rk, rl, rm, rn, rp
         real *16, DIMENSION(5, 6) :: coeff1
         real *16, DIMENSION(5, 6) :: coeff2

         if(set_read.Eq.0) then
         ne = ng/(p+1)
         pp = p+1

         !CARTESIAN GRIDPOINTS
         coeff1(1,1) = -1.0000000000000000D0
         coeff1(1,2) = +1.0000000000000000D0

         coeff1(2,1) = -1.0000000000000000D0
         coeff1(2,2) = +0.0000000000000000D0
         coeff1(2,3) = +1.0000000000000000D0

         coeff1(3,1) = -1.0000000000000000D0
         coeff1(3,2) = -4.47213595499957939D-1
         coeff1(3,3) = +4.47213595499957939D-1
         coeff1(3,4) = +1.0000000000000000D0

         coeff1(4,1) = -1.0000000000000000D0
         coeff1(4,2) = -6.54653670707977144D-1
         coeff1(4,3) = +0.0000000000000000D0
         coeff1(4,4) = +6.54653670707977144D-1
         coeff1(4,5) = +1.0000000000000000D0

         coeff1(5,1) = -1.0000000000000000D0
         coeff1(5,2) = -7.6505532392946469D-1
         !coeff1(5,3) = -2.852315164806e-01
         coeff1(5,3) = -2.85231516480645096D-1
         !coeff1(5,4) = +2.852315164806e-01
         coeff1(5,4) = +2.85231516480645096D-1
         coeff1(5,5) = +7.6505532392946469D-1
         coeff1(5,6) = +1.0000000000000000D0

         !OUTPUT GRIDPOINTS
         coeff2(1,1) = -1.0000000000000000D0
         coeff2(1,2) = +1.0000000000000000D0

         coeff2(2,1) = -1.0000000000000000D0
         coeff2(2,2) = +0.0000000000000000D0
         coeff2(2,3) = +1.0000000000000000D0

         coeff2(3,1) = -1.0000000000000000D0
         coeff2(3,2) = -3.3333333333333333D-1
         coeff2(3,3) = +3.3333333333333333D-1
         coeff2(3,4) = +1.0000000000000000D0

         coeff2(4,1) = -1.0000000000000000D0
         coeff2(4,2) = -5.0000000000000000D-1
         coeff2(4,3) = +0.0000000000000000D0
         coeff2(4,4) = +5.0000000000000000D-1
         coeff2(4,5) = +1.0000000000000000D0

         coeff2(5,1) = -1.0000000000000000D0
         coeff2(5,2) = -6.0000000000000000D-1
         coeff2(5,3) = -2.0000000000000000D-1
         coeff2(5,4) = +2.0000000000000000D-1
         coeff2(5,5) = +6.0000000000000000D-1
         coeff2(5,6) = +1.0000000000000000D0

          if( ng .Gt. MaxG ) then
            write(*,*) ' Error : incorrect size of grid'
            stop
          endif
          write(*,*) ne
          Ni1 = 1
          Nj1 = 1
          Nk1 = 1
          Ni2 = ng
          Nj2 = ng
          Nk2 = ng
          Ni3 = ne
          Nj3 = ne
          Nk3 = ne

          dh = 2.0D0 * Pi_const
          dh = dh/dble(ne)

          do k = Nk1, Nk3, 1
          do j = Nj1, Nj3, 1
          do i = Ni1, Ni3, 1
            do l = 1, pp, 1
            do m = 1, pp, 1
            do n = 1, pp, 1
              
              ri = dble(i)
              rj = dble(j)
              rk = dble(k)
              rl = dble(l)
              rm = dble(m)
              rn = dble(n)
              rp = dble(p)

              indi = (i-1)*pp+l
              indj = (j-1)*pp+m
              indk = (k-1)*pp+n

              !At Cartesian Points
              Xg(indi,indj,indk) = (ri-1)*dh + (coeff2(p,l)+1)/2*dh;
              Yg(indi,indj,indk) = (rj-1)*dh + (coeff2(p,m)+1)/2*dh;
              Zg(indi,indj,indk) = (rk-1)*dh + (coeff2(p,n)+1)/2*dh;

              !At Solution Points
              Xs(indi,indj,indk) = (ri-1)*dh + (coeff1(p,l)+1)/2*dh;
              Ys(indi,indj,indk) = (rj-1)*dh + (coeff1(p,m)+1)/2*dh;
              Zs(indi,indj,indk) = (rk-1)*dh + (coeff1(p,n)+1)/2*dh;

            enddo
            enddo
            enddo
          enddo
          enddo
          enddo
          
          VolBox = dble((Ni2-Ni1)*(Nj2-Nj1)*(Nk2-Nk1))
          Ngrd   = 3

          else

          if( ng .Gt. MaxG ) then
            write(*,*) ' Error : incorrect size of grid'
            stop
          endif
          Ni1 = 1
          Nj1 = 1
          Nk1 = 1
          Ni2 = ng
          Nj2 = ng
          Nk2 = ng
          dh = (2.0 * Pi_const)/dble(ng - 1)
          do k = Nk1, Nk2, 1
          do j = Nj1, Nj2, 1
          do i = Ni1, Ni2, 1
            Xg(i,j,k) = (i-1) * dh
            Yg(i,j,k) = (j-1) * dh
            Zg(i,j,k) = (k-1) * dh
          enddo
          enddo
          enddo

          VolBox = dble((Ni2-Ni1)*(Nj2-Nj1)*(Nk2-Nk1))
          Ngrd   = 3

          endif
        end

