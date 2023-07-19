c***********************************************************************
c    Gauss elimination
c***********************************************************************
c
      program Gauss_elimination
c
      implicit real*8(a-h,o-z)    ! 暗黙の型宣言，倍精度
        ! 暗黙の型宣言を使用しない場合は，「implicit none」として
        ! 変数を定義すること
c
      parameter(N=4)      ! 配列の領域確保のため
c
      dimension A(N,N),X(N),B(N),AL(N,N),AD(N,N),AU(N,N)
c
      open(2,file='rst.txt')
c
      call inp(A,B,N)
c
      call gauss(A,X,B,AL,AD,AU,N)
c
      write(2,*) 'AL'
      do i=1,N
        write(2,'(10E12.4)') AL(i,1:N)
      end do
c
      write(2,*) 'AD'
      do i=1,N
        write(2,'(10E12.4)') AD(i,1:N)
      end do
c
      write(2,*) 'AU'
      do i=1,N
        write(2,'(10E12.4)') AU(i,1:N)
      end do
c
      write(2,*) 'X'
      do i=1,N
        write(2,'(E12.4)') X(i)
      end do
c
      close(2)
c
      stop
      end
c
c
c
c***********************************************************************
c     analytical model input
c***********************************************************************
      subroutine inp(a,b,N)
c
      implicit real*8(a-h,o-z)
c
      dimension a(N,*),b(*)
c
      a(1,1) = 5D0
      a(2,1) = -4D0
      a(3,1) = 1D0
      a(4,1) = 0D0
      a(1,2) = a(2,1)
      a(2,2) = 6D0
      a(3,2) = -4D0
      a(4,2) = 1D0
      a(1,3) = a(3,1)
      a(2,3) = a(3,2)
      a(3,3) = 6D0
      a(4,3) = -4D0
      a(1,4) = a(4,1)
      a(2,4) = a(4,2)
      a(3,4) = a(4,3)
      a(4,4) = 5D0
c
      b(1) = 0D0
      b(2) = 1D0
      b(3) = 0D0
      b(4) = 0D0
c
      return
      end
c
c
c
