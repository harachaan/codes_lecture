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
      dimension D(3,3),TK(2*MAXNODE,2*MAXNODE)
c
      open(2,file='rst.txt')
      open(3,file='inp_hole.txt')
c
      call inp(NNODE,NDOF,NELEM,X,Y,THICK,LNODS,E,VNU)
      close(3)
      call dmat(D,E,VNU)
c
      TK(1:NDOF,1:NDOF) = 0.D0   ! initialize
c
      do lelem=1,NELEM
        call bmat(lelem,B,X,Y,LNODS)
        call mergekmat(TK,lelem,LNODS,B,X,Y,D,THICK)
      end do
c
      write(2,*) 'TK(i,i)'
      do i=1,NDOF
        write(2,'(E15.7)') TK(i,i)
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
      subroutine inp(nnode,ndof,nelem,x,y,thick,lnods,e,vnu)
c
      implicit real*8(a-h,o-z)
c
      dimension x(*),y(*),lnods(3,*)
c
c---- �@�B�I����
      e = 200D3
      vnu = 0.4D0
c
c---- �����̒�`
      thick = 1.D0
c
c---- �C���v�b�g�t�@�C�����
      nnode = 30      ! �ߓ_��
      nelem = 40      ! �v�f��
c
c---- �ߓ_���W
      do i=1,nnode
        read(3,*) inode,x(inode),y(inode)
c        write(2,*) x(inode),y(inode)
      end do
c
      ndof = 2 * nnode   ! �����R�x���i���m�ψʐ����̐��j
c
c---- �v�f���\������ߓ_�̃��X�g�i�K�������v���j
      do l=1,nelem
        read(3,*) lelem,lnods(1:3,lelem)
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
