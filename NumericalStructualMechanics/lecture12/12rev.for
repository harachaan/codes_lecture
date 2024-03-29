c***********************************************************************
c     Finite Element Analysis
c     Element type: 4-node quadrilateral elements (2-dimensional)
c***********************************************************************
c
      program FEM_quad_element
c
      implicit real*8(a-h,o-z)    ! 暗黙の型宣言，倍精度
        ! 暗黙の型宣言を使用しない場合は，「implicit none」として
        ! 変数を定義すること
c
      include 'param.inc'
c
      dimension X(MAXNODE),Y(MAXNODE),LNODS(4,MAXELEM),B(3,8,4)
      dimension D(3,3),TK(2*MAXNODE,2*MAXNODE),U(2*MAXNODE),DJ(4)
      dimension P(2*MAXNODE),I_BC_GIVEN(MAXBC),V_BC_GIVEN(MAXBC)
      dimension STRAIN(3,4,MAXELEM),STRESS(3,4,MAXELEM)
      dimension R1(4),R2(4),DNDR(2,4,4),w1(4),w2(4)
      dimension sgm(3,MAXNODE),nsum(MAXNODE)
c
      open(2,file='rst.txt')
      open(3,file='inp_concentrated-f_quadrilateral.txt')
c
      call inp(NNODE,NDOF,NELEM,X,Y,THICK,LNODS,E,VNU,
     &         P,N_BC_GIVEN,I_BC_GIVEN,V_BC_GIVEN)
      close(3)
      call mapping(r1,r2,dndr,w1,w2)
      call dmat(D,E,VNU)
c
      TK(1:NDOF,1:NDOF) = 0.D0   ! initialize
c
      do lelem=1,NELEM
        call bmat(lelem,B,DJ,X,Y,LNODS,DNDR)
        call mergekmat(TK,lelem,LNODS,B,DJ,D,THICK,W1,W2)
      end do
c
      call bound(TK,P,NDOF,N_BC_GIVEN,I_BC_GIVEN,V_BC_GIVEN)
c
      call gauss2(TK,U,P,NDOF)
c
      do lelem=1,NELEM
        call bmat(lelem,B,DJ,X,Y,LNODS,DNDR)
        call calc_stress(lelem,STRAIN,STRESS,B,D,U,LNODS)
      end do
c
      write(2,*) 'STRESS'
      do lelem=1,NELEM
        write(2,'(I5,A8,3E15.7)') lelem,'integ1',STRESS(1:3,1,lelem)
        write(2,'(I5,A8,3E15.7)') lelem,'integ2',STRESS(1:3,2,lelem)
        write(2,'(I5,A8,3E15.7)') lelem,'integ3',STRESS(1:3,3,lelem)
        write(2,'(I5,A8,3E15.7)') lelem,'integ4',STRESS(1:3,4,lelem)
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
      subroutine inp(nnode,ndof,nelem,x,y,thick,lnods,e,vnu,
     &                p,n_bc_given,i_bc_given,v_bc_given)
c
      implicit real*8(a-h,o-z)
c
      dimension x(*),y(*),lnods(4,*),p(*),i_bc_given(*),v_bc_given(*)
c
c---- 機械的特性
      e = 200D3
      vnu = 0.4D0
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
        read(3,*) lnods(1:4,lelem)
      end do
c
c      do lelem=1,nelem
c        write(2,*) lnods(1:4,lelem)
c      end do
c
c---- 境界条件の設定
c---- 変位境界条件（x方向）
      num = 0     ! 境界条件の番号，インクリメントする
      idof = 1    ! x方向
      do inode=1,nnode
        if((x(inode) .eq. 0D0) .and. (y(inode) .eq. 0D0)) then
          num = num + 1
          i_bc_given(num) = inode     ! 節点番号
          v_bc_given(num) = 0.D0      ! 強制変位量
          i_bc_given(num) = (i_bc_given(num) - 1) * 2 + idof
        end if
      end do
c
c---- 変位境界条件（y方向）
      idof = 2    ! y方向
      do inode=1,nnode
        if(y(inode) .eq. 0D0) then
          num = num + 1
          i_bc_given(num) = inode     ! 節点番号
          v_bc_given(num) = 0.D0      ! 強制変位量
          i_bc_given(num) = (i_bc_given(num) - 1) * 2 + idof
        end if
      end do
c
      n_bc_given = num
c
c---- 力学的境界条件
      p(1:ndof) = 0.d0    ! 初期化（力の作用しない自由表面では，節点力はゼロ）
c
      idof = 2                  ! y方向
      do inode=1,nnode
        if((x(inode) .eq. 0D0) .and. (y(inode) .eq. 50D0)) then
          i_load_given = (inode - 1) * 2 + idof   ! 自由度番号
          p(i_load_given) = - 10D0 * thick   ! 集中荷重10N
          exit
        end if
      end do
c
      return
      end
c
c
c
