c***********************************************************************
c     Finite Element Analysis
c     Element type: 3-node triangle elements (2-dimensional)
c***********************************************************************
c
      program FEM_triangle_element
c
      implicit real*8(a-h,o-z)    ! 暗黙の型宣言，倍精度
        ! 暗黙の型宣言を使用しない場合は，「implicit none」として
        ! 変数を定義すること
c
      include 'param.inc'
c
      dimension X(MAXNODE),Y(MAXNODE),LNODS(3,MAXELEM),B(3,6)
c
      open(2,file='rst.txt')
      open(3,file='inp_concentrated-f_triangle.txt')
c
      call inp(NNODE,NDOF,NELEM,X,Y,THICK,LNODS)
      close(3)
c
      do lelem=1,NELEM
        call bmat(lelem,B,X,Y,LNODS)
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
      subroutine inp(nnode,ndof,nelem,x,y,thick,lnods)
c
      implicit real*8(a-h,o-z)
c
      dimension x(*),y(*),lnods(3,*)
c
c---- 厚さの定義
      thick = 1.D0
c
c---- インプットファイルより節点数，要素数
      read(3,*) nnode,nelem
c
c---- 節点座標
      do inode=1,nnode
        read(3,*) x(inode),y(inode)
c        write(2,*) x(inode),y(inode)
      end do
c
      ndof = 2 * nnode   ! 総自由度数（未知変位成分の数）
c
c---- 要素を構成する節点のリスト（必ず反時計回り）
      do lelem=1,nelem
        read(3,*) lnods(1:3,lelem)
      end do
c
c      do lelem=1,nelem
c        write(2,*) lnods(1:3,lelem)
c      end do
c
      return
      end
c
c
c
