! **********************************************************************
!    Finite Element Analysis
!    Element type: 4-node quadrilateral elements (2-dimensional)
! **********************************************************************
! 
program FEM_quad_element

    implicit real*8(a-h,o-z)    ! 暗黙の型宣言，倍精度
    ! 暗黙の型宣言を使用しない場合は，「implicit none」として
    ! 変数を定義すること

    include '../assignment3/parameter/param.inc'

    dimension X(MAXNODE),Y(MAXNODE),LNODS(4,MAXELEM),B(3,8,4)
    dimension D(3,3),TK(2*MAXNODE,2*MAXNODE),U(2*MAXNODE),DJ(4)
    dimension P(2*MAXNODE),I_BC_GIVEN(MAXBC),V_BC_GIVEN(MAXBC)
    dimension STRAIN(3,4,MAXELEM),STRESS(3,4,MAXELEM)
    dimension R1(4),R2(4),DNDR(2,4,4),w1(4),w2(4)
    ! dimension sgm(3,MAXNODE),nsum(MAXNODE)

    open(2,file='rst.txt')
    open(3,file='../assignment3/parameter/inp_concentrated-f_quadrilateral.txt')

    call inp(NNODE, NDOF, NELEM, X, Y, THICK, LNODS, E, VNU, P, N_BC_GIVEN, I_BC_GIVEN, V_BC_GIVEN)
    close(3)
    call mapping(r1, r2, dndr, w1, w2)
    call dmat(D, E, VNU)

    TK(1:NDOF, 1:NDOF) = 0.D0   ! initialize

    do lelem=1, NELEM
        call bmat(lelem, B, DJ, X, Y, LNODS, DNDR)
        call mergekmat(TK, lelem, LNODS, B, DJ, D, THICK, W1, W2, MAXNODE)
    end do

    call bound(TK, P, NDOF, N_BC_GIVEN, I_BC_GIVEN, V_BC_GIVEN, MAXNODE)

    call gauss2(TK, U, P, NDOF, MAXNODE)


    write(2, *) 'Bmat'
    do lelem=1, NELEM
        call bmat(lelem, B, DJ, X, Y, LNODS, DNDR)
        call calc_stress(lelem, STRAIN, STRESS, B, D, U, LNODS)


        if (lelem == 1) then         
        do integ=1,4             
            write(2,*) 'integ=',integ             
            do i=1,3                 
                do j=1,8                     
                    if (j .eq. 8) then                         
                        write(2, '(E12.4)') B(i,j,integ)                     
                    else                         
                        write(2, '(E12.4)',advance='no') B(i,j,integ)                     
                    end if                 
                end do             
            end do         
        end do     
        end if
    end do



    write(2, *) 'STRESS'
    do lelem = 1, NELEM
        write(2, '(I5, A8, 3E15.7)') lelem, 'integ1', STRESS(1:3, 1, lelem)
        write(2, '(I5, A8, 3E15.7)') lelem, 'integ2', STRESS(1:3, 2, lelem)
        write(2, '(I5, A8, 3E15.7)') lelem, 'integ3', STRESS(1:3, 3, lelem)
        write(2, '(I5, A8, 3E15.7)') lelem, 'integ4', STRESS(1:3, 4, lelem)
    end do

    write(2, *) 'Paste to Excel -> lelem=41 to 50'
    ! choose lelem=41 to 50 or lelem=51 to 60. Both are same meaning.
    write(2, '(I5, A8, 3E15.7)') 41, 'integ2', STRESS(1:3, 2, 41)        ! 始点
    do lelem = 41, 50
        write(2, '(I5, A8, 3E15.7)') lelem, 'integ3', STRESS(1:3, 3, lelem)
    end do

    write(2, *) 'Paste to Excel -> lelem=51 to 60'
    write(2, '(I5, A8, 3E15.7)') 51, 'integ1', STRESS(1:3, 1, 51)        ! 始点
    do lelem = 51, 60
        write(2, '(I5, A8, 3E15.7)') lelem, 'integ4', STRESS(1:3, 4, lelem)
    end do

    close(2)

    stop
    end



! **********************************************************************
!     analytical model input
! **********************************************************************
subroutine inp(nnode,ndof,nelem,x,y,thick,lnods,e,vnu,p,n_bc_given,i_bc_given,v_bc_given)

    implicit real*8(a-h,o-z)

    dimension x(*),y(*),lnods(4,*),p(*),i_bc_given(*),v_bc_given(*)

!---- 機械的特性
    e = 200D3
    vnu = 0.4D0

!---- 厚さの定義
    thick = 1.D0

!---- インプットファイルより節点数，要素数
    read(3,*) nnode,nelem

!---- 節点座標
    do inode=1,nnode
        read(3,*) x(inode),y(inode)
        ! write(2,*) x(inode),y(inode)
    end do

    ndof = 2 * nnode   ! 総自由度数（未知変位成分の数）

!---- 要素を構成する節点のリスト（必ず反時計回り）
    do lelem=1,nelem
        read(3,*) lnods(1:4,lelem)
    end do

!---- 境界条件の設定
!---- 変位境界条件（x方向）
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

!---- 変位境界条件（y方向）
    idof = 2    ! y方向
    do inode=1,nnode
    if(y(inode) .eq. 0D0) then
        num = num + 1
        i_bc_given(num) = inode     ! 節点番号
        v_bc_given(num) = 0.D0      ! 強制変位量
        i_bc_given(num) = (i_bc_given(num) - 1) * 2 + idof
    end if
    end do

    n_bc_given = num

!---- 力学的境界条件
    p(1:ndof) = 0.D0    ! 初期化（力の作用しない自由表面では，節点力はゼロ）

    idof = 2                  ! y方向
    do inode=1,nnode
    if((x(inode) .eq. 0D0) .and. (y(inode) .eq. 50D0)) then
        i_load_given = (inode - 1) * 2 + idof   ! 自由度番号
        p(i_load_given) = - 10D0 * thick   ! 集中荷重10N
        exit
    end if
    end do

    return
end


subroutine bmat(lelem, B, dj, x, y, lnods, dndr)

    implicit real*8(a-h,o-z)
    dimension x(*), y(*), lnods(4, *), B(3, 8, 4)
    dimension dj(4), xe(2, 4), dxdr(2, 2), dndx(2, 4), dndr(2, 4, 4), temp(2, 2)

    b(1:3, 1:8, 1:4) = 0.D0

    do n = 1, 4
        xe(1, n) = x(lnods(n, lelem))
        xe(2, n) = y(lnods(n, lelem))
    end do

    do integ = 1, 4
        do i = 1, 2
            do j = 1, 2
                dxdr(j, i) = 0.D0
                do k = 1, 4
                    dxdr(j, i) = dxdr(j, i) + dndr(j, k, integ) * xe(i, k)
                end do
            end do
        end do

        dj(integ) = dxdr(1, 1) * dxdr(2, 2) - dxdr(1, 2) * dxdr(2, 1)
        temp = dxdr
        dxdr(1, 1) = temp(2, 2)/dj(integ)
        dxdr(1, 2) = -temp(1, 2)/dj(integ)
        dxdr(2, 1) = -temp(2, 1)/dj(integ)
        dxdr(2, 2) = temp(1, 1)/dj(integ)


        do n = 1, 4
            do i = 1, 2
                dndx(i, n) = 0.D0
                do k = 1, 2
                    dndx(i, n) = dndx(i, n) + dxdr(i, k) * dndr(k, n, integ)
                end do
            end do
        end do

        do n = 1, 4
            i1 = (n-1) * 2 + 1
            i2 = (n-1) * 2 + 2
            B(1, i1, integ) = dndx(1, n)
            B(2, i2, integ) = dndx(2, n)
            B(3, i1, integ) = dndx(2, n)
            B(3, i2, integ) = dndx(1, n)
        end do
    end do

    return
end


subroutine dmat(D, E, vnu)
    implicit real*8(a-h,o-z)
    dimension D(3, 3)

    ! Dmatrixを初期化
    D(1:3, 1:3) = 0.D0
    ! Dmatrixを計算
    ! 平面ひずみ -> eq.(6)
    D(1, 1) = 1.0D0 - vnu
    D(1, 2) = vnu
    D(2, 1) = D(1, 2)
    D(2, 2) = 1.0D0 - vnu
    D(3, 3) = (1.0D0 - 2.0D0*vnu) / 2.0D0
    D = D * E / ((1.0D0+vnu) * (1.0D0-2.0D0*vnu))
    return
end


subroutine mergekmat(tk, lelem, lnods, b, dj, d, thick, w1, w2, MAXNODE)
    implicit real*8(a-h,o-z)

    dimension ip(8), lnods(4,*)
    dimension tk(2*MAXNODE, 2*MAXNODE), b(3, 8, 4), d(3, 3), db(3, 8), dk(8, 8), dj(4), dtdb(8, 8), w1(*), w2(*)

    dk(1:8, 1:8) = 0.D0
    ! ip(1:6) = 0.D0      ! initialize

    do inode = 1, 4
        node_no = lnods(inode, lelem)
        do idof = 1, 2
            ip((inode-1)*2+idof) = (node_no-1) * 2 + idof
        end do
    end do

    do integ = 1, 4
        do i = 1, 3
            do j = 1, 8
                db(i, j) = 0.D0     ! initialize
                do k = 1, 3
                    db(i, j) = db(i, j) + d(i, k) * b(k, j, integ)
                end do
            end do
        end do

        do i = 1, 8
            do j = 1, 8
                dtdb(i, j) = 0.D0     ! initialize
                do k = 1, 3
                    dtdb(i, j) = dtdb(i, j) + b(k, i, integ) * db(k, j) * dj(integ)
                end do
            end do
        end do

        do i = 1, 8
            do j = 1, 8
                ! dk(i, j) = 0.D0     ! initialize
                dk(i, j) = dk(i, j) + dtdb(i, j) * thick * w1(integ) * w2(integ)
            end do
        end do
    end do


    do i = 1, 8
        do j = 1, 8
            tk(ip(i), ip(j)) = tk(ip(i), ip(j)) + dk(i, j)
        end do
    end do
    
    return 
end


subroutine bound(TK, P, NDOF, N_BC_GIVEN, I_BC_GIVEN, V_BC_GIVEN, MAXNODE)
    implicit real*8(a-h,o-z)
    ! 境界条件を設定する関数
    ! N_BC_GIVEN: 変位境界条件の数
    ! I_BC_GIVEN: 変位境界条件を設定する自由度番号
    ! V_BC_GIVEN: 変位境界条件の既知変位の値

    integer MAXNODE
    dimension TK(2*MAXNODE,2*MAXNODE), P(*), I_BC_GIVEN(*), V_BC_GIVEN(*)
    integer NDOF, N_BC_GIVEN


    ! NDOF 総自由度数
    do i = 1, N_BC_GIVEN
        j = I_BC_GIVEN(i)
        do k = 1, NDOF
            P(k) = P(k) - V_BC_GIVEN(i) * TK(k, j)
        end do
    end do

    do i = 1, N_BC_GIVEN
        j = I_BC_GIVEN(i)
        P(j) = V_BC_GIVEN(i)
    end do

    do i = 1, N_BC_GIVEN
        j = I_BC_GIVEN(i)
        do k = 1, NDOF
            TK(k, j) = 0.D0
            TK(j, k) = 0.D0
        end do
        TK(j, j) = 1.D0
    end do
    return
end


subroutine gauss2(A, y, B, N, MAXNODE)
    ! 配列を合理化したgauss
    implicit real*8(a-h, o-z)
    integer N
    double precision A(2*MAXNODE,2*MAXNODE), B(2*MAXNODE), y(2*MAXNODE)

    A(1, 2) = A(1, 2) / A(1, 1)
    A(2, 1) = A(2, 1) / A(1, 1)
    A(2, 2) = A(2, 2) - A(2, 1) * A(1, 1) * A(1, 2)

    do j = 3, N
        A(1, j) = A(1, j) / A(1, 1)
        A(j, 1) = A(j, 1) / A(1, 1)

        do i = 2, j-1
            ! 上三角の2行目から対角の前まで
            sum_U = 0.D0 ! initialize
            sum_L = 0.D0
            do k = 1, i-1
                sum_U = sum_U + A(i, k) * A(k, k) * A(k, j)
                sum_L = sum_L + A(j, k) * A(k, k) * A(k, i)
            end do
                
            A(i, j) = (A(i, j)-sum_U) / A(i, i)
            A(j, i) = (A(j, i)-sum_L) / A(i, i)
        end do

        sum_D = 0.D0
        do k = 1, j-1
            sum_D = sum_D + A(j, k) * A(k, k) * A(k, j)
        end do
        A(j, j) = A(j, j) - sum_D
    end do

    ! 前進代入
    y = 0.D0
    y(1) = B(1)
    do i = 2, N
        temp = 0.D0        ! intialize
        do j = 1, i-1
            temp = temp + A(i, j) * y(j)
        end do
        y(i) = B(i) - temp
    end do

    ! 対角項
    do i = 1, N
        y(i) = y(i) / A(i, i)
    end do

    ! 後進消去
    do j = N, 2, -1
        do k = 1, j-1
            y(k) = y(k) - A(k, j) * y(j)
        end do
    end do
    return
end


subroutine calc_stress(lelem, strain, stress, b, d, u, lnods)
    implicit real*8(a-h,o-z)
    dimension strain(3, *), stress(3, *), b(3, *), d(3, *), u(*), lnods(3, *)

    strain(1:3, lelem) = 0.D0
    do i = 1, 3
        do j = 1, 3
            ! 節点番号 -> lelem=1のとき、lnods=(1, 12, 2) のいずれか
            inode = lnods(j, lelem)
            ! eq. (3.15): strain = Bmat * 節点変位ベクトル
            strain(i, lelem) = strain(i, lelem) + b(i, 2*j-1)*u(2*inode-1) + b(i, 2*j)*u(2*inode)
        end do
    end do

    stress(1:3, lelem) = 0.D0
    do i = 1, 3
        do j = 1, 3
            ! eq. (4.14): stress = Dmat * strain
            stress(i, lelem) = stress(i, lelem) + d(i, j) * strain(j, lelem)
        end do
    end do
    return    
end


subroutine mapping(r1, r2, dndr, w1, w2)
    implicit real*8(a-h,o-z)
    dimension r1(4), r2(4), dndr(2, 4, 4), w1(4), w2(4)

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


    do integ = 1, 4
        dndr(1, 1, integ) = -0.25D0 * (1.D0 - r2(integ))
        dndr(2, 1, integ) = -0.25D0 * (1.D0 - r1(integ))
        dndr(1, 2, integ) =  0.25D0 * (1.D0 - r2(integ))
        dndr(2, 2, integ) = -0.25D0 * (1.D0 + r1(integ))
        dndr(1, 3, integ) =  0.25D0 * (1.D0 + r2(integ))
        dndr(2, 3, integ) =  0.25D0 * (1.D0 + r1(integ))
        dndr(1, 4, integ) = -0.25D0 * (1.D0 + r2(integ))
        dndr(2, 4, integ) =  0.25D0 * (1.D0 - r1(integ))
    end do



    return
end