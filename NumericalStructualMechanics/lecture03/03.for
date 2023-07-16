c***********************************************************************
c     Finite Element Analysis
c     Element type: 3-node triangle elements (2-dimensional)
c***********************************************************************
c
      program FEM_triangle_element
c
      implicit real*8(a-h,o-z)    ! �Öق̌^�錾�C�{���x
        ! �Öق̌^�錾���g�p���Ȃ��ꍇ�́C�uimplicit none�v�Ƃ���
        ! �ϐ����`���邱��
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
c---- �����̒�`
      thick = 1.D0
c
c---- �C���v�b�g�t�@�C�����ߓ_���C�v�f��
      read(3,*) nnode,nelem
c
c---- �ߓ_���W
      do inode=1,nnode
        read(3,*) x(inode),y(inode)
c        write(2,*) x(inode),y(inode)
      end do
c
      ndof = 2 * nnode   ! �����R�x���i���m�ψʐ����̐��j
c
c---- �v�f���\������ߓ_�̃��X�g�i�K�������v���j
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
