program homework8


implicit real*8(a-h,o-z)

include 'param.inc'

dimension X(MAXNODE), Y(MAXNODE), LNODS(3,MAXELEM), B(3,8,4)
  !MAXNODE = 1000
dimension D(3,3),TK(2*MAXNODE,2*MAXNODE),U(2*MAXNODE),DJ(4)
dimension P(2*MAXNODE),I_BC_GIVEN(MAXBC),V_BC_GIVEN(MAXBC)
  !P>>右辺ベクトル,I_BC~>>変異境界条件を設定する自由度,V_BC~>>変異境界条件の既知変異の値
dimension STRAIN(3,4,MAXELEM),STRESS(3,4,MAXELEM)
dimension R1(4),R2(4),DNDR(2,4,4),w1(4),w2(4)
dimension sgm(3,MAXNODE),nsum(MAXNODE)

open(3,file='inp_concentrated-f_quadrilateral.txt')
open(4,file='homework8_result_梶川拓海.txt')

call inp(NNODE,NDOF,NELEM,X,Y,THICK,LNODS,E,VNU,P,N_BC_GIVEN,I_BC_GIVEN,V_BC_GIVEN)
! NNODE=121, NDOF=242, NELEM=200, 
call mapping(r1,r2,dndr,w1,w2)
call dmat(D,E,VNU)

!全体剛性マトリックス
TK(1:NDOF,1*NDOF) = 0.D0

do lelem = 1,NELEM
    call bmat(lelem,B,DJ,X,Y,LNODS,DNDR)
    call mergekmat(TK,lelem,LNODS,B,DJ,D,thick,W1,W2)
end do

call bound(TK,P,NDOF,N_BC_GIVEN,I_BC_GIVEN,V_BC_GIVEN)

call gauss2(TK,U,P,NDOF)

do lelem = 1,NELEM
    call bmat(lelem,B,DJ,X,Y,LNODS,DNDR)
    call calc_stress(lelem,STRAIN,STRESS,B,D,U,LNODS)
end do

write(4,*) 'STRESS'
do lelem=1,NELEM
    write(4,'(I5,A8,3E15.7)') lelem,'integ1',STRESS(1:3,1,lelem)
    write(4,'(I5,A8,3E15.7)') lelem,'integ2',STRESS(1:3,2,lelem)
    write(4,'(I5,A8,3E15.7)') lelem,'integ3',STRESS(1:3,3,lelem)
    write(4,'(I5,A8,3E15.7)') lelem,'integ4',STRESS(1:3,4,lelem)
end do

close(3)
close(4)

stop
end



subroutine inp(nnode,ndof,nelem,x,y,thick,lnods,e,vnu,p,n_bc_given,i_bc_given,v_bc_given)

implicit real*8(a-h,o-z)
dimension x(*),y(*),lnods(4,*),p(*),i_bc_given(*),v_bc_given(*)

!機械的特性
e = 200D3
vnu = 0.4D0

thick = 1.D0

read(3,*) nnode,nelem !節点数，要素数

do inode = 1,nnode
    read(3,*) x(inode),y(inode)
    !write(2,*) x(inode),y(inode)
end do

ndof = 2 * nnode !総自由度数

do lelem = 1,nelem
    read(3,*) lnods(1:4,lelem)
end do

!変異境界条件(x方向)
num = 0 !境界条件の番号
idof = 1 !x方向
do inode = 1,nnode
    if((x(inode) == 0.D0) .and. (y(inode) == 0.D0))then
        num = num +1
        i_bc_given(num) = inode!節点番号
        v_bc_given(num) = 0.D0 !強制変位量
        i_bc_given(num) = (i_bc_given(num) - 1) * 2 + idof
    end if
end do

!変異境界条件(y方向)
idof = 2 !y方向
do inode = 1,nnode
    if(y(inode) == 0.D0) then
        num = num +1
        i_bc_given(num) = inode
        v_bc_given(num) = 0.D0 !強制変位量
        i_bc_given(num) = (i_bc_given(num) - 1) * 2 + idof
    end if
end do

n_bc_given = num

!力学的境界条件
p(1:ndof) = 0.D0 !初期化

idof = 2 !y方向
do inode = 1, nnode
    if((x(inode) == 0.D0) .and. (y(inode) == 50.D0)) then
        i_load_given = (inode-1) * 2 + idof !自由度番号
        p(i_load_given) = -10.D0 * thick !集中荷重
        exit
    end if
end do

return 
end



subroutine mapping(r1,r2,dndr,w1,w2)

implicit real*8(a-h,o-z)
dimension r1(4),r2(4),dndr(2,4,4),w1(4),w2(4)

r1(1) = -1.D0
r1(2) = 1.D0
r1(3) = 1.D0
r1(4) = -1.D0
r1(1:4) = r1(1:4)/dsqrt(3.D0)

w1(1:4) = 1.D0

r2(1) = -1.D0
r2(2) = -1.D0
r2(3) = 1.D0
r2(4) = 1.D0
r2(1:4) = r2(1:4)/dsqrt(3.D0)

w2(1:4) = 1.D0

do i = 1,4
    dndr(1,1,i) = -0.25D0*(1.D0-r2(i))
    dndr(2,1,i) = -0.25D0*(1.D0-r1(i))
    dndr(1,2,i) = 0.25D0*(1.D0-r2(i))
    dndr(2,2,i) = -0.25D0*(1.D0+r1(i))
    dndr(1,3,i) = 0.25D0*(1.D0+r2(i))
    dndr(2,3,i) = 0.25D0*(1.D0+r1(i))
    dndr(1,4,i) = -0.25D0*(1.D0+r2(i))
    dndr(2,4,i) = 0.25D0*(1.D0-r1(i))
end do

return
end



subroutine bmat(lelem,b,dj,x,y,LNODS,dndr)

implicit real*8(a-h,o-z)
dimension xe(2,4),LNODS(4,*),x(*),y(*),b(3,8,4)
dimension dndr(2,4,4),dj(4),dxdr(2,2),dndx(2,4),dxdr2(2,2)


xe(1:2,1:4) = 0.D0
b(1:3,1:8,1:4) = 0.D0
dxdr2(1:2,1:2) = 0.D0

do m = 1,4
    inode = LNODS(m,LELEM)
    xe(1,m) = X(inode)
    xe(2,m) = Y(inode)
end do

do integ = 1,4
    do i = 1,2
        do j = 1,2
            dxdr(j,i) = 0.D0
            do k = 1,4
                dxdr(j,i) = dxdr(j,i) + dndr(j,k,integ) * xe(i,k)
            end do
        end do
    end do

    dj(integ) = dxdr(1,1)*dxdr(2,2)-dxdr(1,2)*dxdr(2,1) !det[dxdr]
    dxdr2(1,1) = dxdr(2,2)/dj(integ)   !dxdrの逆行列
    dxdr2(1,2) = -dxdr(1,2)/dj(integ)
    dxdr2(2,1) = -dxdr(2,1)/dj(integ)
    dxdr2(2,2) = dxdr(1,1)/dj(integ)

    do i = 1,4
        do j = 1,2
            dndx(j,i) = 0.D0
            do k = 1,2
                dndx(j,i) = dndx(j,i) + dxdr2(j,k) * dndr(k,i,integ)
            end do
        end do
    end do

    do i = 1,4
        m = (i-1)*2+1
        n = (i-1)*2+2
        b(1,m,integ) = dndx(1,i)
        b(2,n,integ) = dndx(2,i)
        b(3,m,integ) = dndx(2,i)
        b(3,n,integ) = dndx(1,i)
    end do
end do

if(lelem==1) then
    do integ=1,4
        write(4,*) "integ=", integ
        do i=1,3
            do j=1,8
                if(j==8) then
                    write(4,"(E12.4)") B(i,j,integ)
                else
                    write(4,"(E12.4)",advance='no') B(i,j,integ)
                end if
            end do
        end do
    end do
end if

return
end


subroutine dmat(d,e,vnu)

implicit real*8(a-h,o-z)
dimension d(3,3)

Econst = e/((1.0D0+vnu)*(1.0D0-2.0D0*vnu))

d(1:3,1:3) = 0.D0

do i = 1,2
    do j = 1,2
        if (i .eq. j) then
            d(i,j) = 1.D0-vnu
        else
            d(i,j) = vnu
        end if
    end do
end do

d(3,3) = (1.0D0-2.D0*vnu)/2.D0

d = Econst*d

return
end


subroutine mergekmat(TK,lelem,lnods,B,dj,d,thick,w1,w2)
implicit real*8(a-h,o-z)
dimension TK(2*1000,2*1000),B(3,8,4),d(3,3),db(3,8,4),dk(8,8),lnods(4,*),ip(8)
dimension w1(4),w2(4),dj(4)

dk(1:8,1:8) = 0.D0

do inode = 1,4
    node_n = lnods(inode,lelem)
    do idof = 1,2
        ip((inode-1)*2+idof) = (node_n-1)*2 + idof
    end do
end do

do integ = 1,4
    do i = 1,3
        do j = 1,8
            db(i,j,integ) = 0.D0
            do k = 1,3
                db(i,j,integ) = db(i,j,integ) + d(i,k) * b(k,j,integ)
            end do
        end do
    end do
    
    do i = 1,8
        do j = 1,8
            dk(i,j) = 0.D0
            do k = 1,3
                dk(i,j) = dk(i,j) + db(k,j,integ) * B(k,i,integ)
            end do
        end do
    end do

    do i = 1,8
        do j = 1,8
            dk(i,j) = w1(integ)*w2(integ)*dk(i,j)*thick*dj(integ)
        end do
    end do

    do i = 1,8
        do j = 1,8
            tk(ip(i),ip(j)) = tk(ip(i),ip(j)) + dk(i,j)
        end do
    end do
end do

return
end


subroutine bound(tk,p,ndof,n_bc_given,i_bc_given,v_bc_given)
implicit real*8(a-h,o-z)
dimension tk(2*1000,2*1000),p(*),i_bc_given(*),v_bc_given(*)

do i = 1,n_bc_given
    j = i_bc_given(i)
    do k = 1,ndof !ndof>>総自由度数
        p(k) = p(k) - v_bc_given(i) * tk(k,j)
    end do
end do

do i = 1,n_bc_given
    j = i_bc_given(i)
    p(j) = v_bc_given(i)
end do

do i = 1,n_bc_given !係数行列いじる
    j = i_bc_given(i)
    do k = 1,ndof
        tk(k,j) = 0
        tk(j,K) = 0
    end do
    tk(j,j) = 1.D0
end do

return
end

subroutine gauss2(tk,u,p,ndof)
implicit real*8(a-h,o-z)
dimension tk(2*1000,2*1000),u(*),p(*),temp(ndof)
dimension sum_U(ndof,ndof),sum_D(ndof,ndof),sum_L(ndof,ndof)
!初期化
sum_U(1:ndof,1:ndof) = 0D0
sum_D(1:ndof,1:ndof) = 0D0
sum_L(1:ndof,1:ndof) = 0D0
temp(1:ndof) = 0D0

tk(1,2) = tk(1,2)/tk(1,1)
tk(2,1) = tk(2,1)/tk(1,1)
tk(2,2) = tk(2,2)-tk(2,1)*tk(1,1)*tk(1,2)

do j = 3,ndof
    tk(1,j) = tk(1,j)/tk(1,1)
    tk(j,1) = tk(j,1)/tk(1,1)
    
    do i = 2,j-1
        do k = 1,i-1
            sum_U(i,j) = sum_U(i,j) + tk(i,k)*tk(k,k)*tk(k,j)
            sum_L(i,j) = sum_L(i,j) + tk(j,k)*tk(k,k)*tk(k,i)
        end do
        
        tk(i,j) = (tk(i,j)-sum_U(i,j))/tk(i,i)
        tk(j,i) = (tk(j,i)-sum_L(i,j))/tk(i,i)
    end do

    do k = 1,j-1
        sum_D(i,j) = sum_D(i,j) + tk(j,k)*tk(k,k)*tk(k,j)
    end do

    tk (j,j) = tk(j,j) - sum_D(i,j)
end do

do i = 2,ndof
    do j = 1,i-1
        temp(i) = temp(i) + tk(i,j)*p(j)
    end do
    p(i) = p(i) - temp(i)
end do

do i = 1,ndof
    p(i) = p(i)/tk(i,i)
end do

u(ndof) = p(ndof)
!後退消去
do i = ndof,2,-1
    do j = 1,i-1
        p(j) = p(j) - tk(j,i)*u(i)
    end do
    u(i-1) = p(i-1)
end do

return
end


subroutine calc_stress(lelem,strain,stress,b,d,u,lnods)
implicit real*8(a-h,o-z)
dimension b(3,8,4),d(3,3),u(*),strain(3,4,*),stress(3,4,*),lnods(4,*)
dimension u_e(8)
!↑節点変位ベクトル

do i = 1,4
    u_e(2*i-1) = u(2*lnods(i,lelem)-1)
    u_e(2*i) = u(2*lnods(i,lelem))
end do

do integ = 1,4
    do i = 1,3
        strain(i,integ,lelem) = 0.D0
        do j = 1,8
            strain(i,integ,lelem) = strain(i,integ,lelem) + b(i,j,integ) * u_e(j)
        end do
    end do
end do

do integ = 1,4
    do i = 1,3
        stress(i,integ,lelem) = 0.D0
        do j = 1,3
            stress(i,integ,lelem) = stress(i,integ,lelem) + d(i,j) * strain(j,integ,lelem)
        end do
    end do
end do

return
end


