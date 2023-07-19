! ---------------------------------------------------------------------
!   Finite Element Analysis
!   Element type: 4-node triangle elements (2-dimensional)
!   4接点の場合の，節点変位ベクトルと各要素の応力を出力する．
!   x=0のy軸に沿ったy方向垂直応力も出力
!
!   Directory Structure
!   - inputFile
!       - inp_concentrated-f_quadrilateral.txt
!       - inp_concentrated-f_triangle.txt
!       - param.inc
!   - lecture12 (curdir)
!       - finalReport.f90
!       - a.out
!       - rst.txt
!       - element_coordinates.txt
!       - kmat.txt
!       - Kmat_diag.txt
!       - dmat.txt
!       - bmat.txt
!       - bound.txt
! ---------------------------------------------------------------------

program FEM_quad_element
    implicit none

    integer, parameter :: DIM = 4

    integer :: MAXNODE, MAXELEM, MAXBC
    include '../inputFile/param.inc' ! MAXNODE, MAXELEM, MAXBCをinclude

    integer :: NNODE, NDOF, NELEM, LNODS(DIM, MAXELEM), lelem, i, N_BC_GIVEN, inode
        ! NNODE: 節点数，NDOF: 総自由度数，NELEM: 要素数, LNODS: 要素を構成する接点のリスト, lelem: 対象要素の番号
    double precision :: X(MAXNODE), Y(MAXNODE), B(3, 8, 4), D(3, 3), TK(2*MAXNODE, 2*MAXNODE), E, VNU, THICK, P(2*MAXNODE), I_BC_GIVEN(MAXBC), V_BC_GIVEN(MAXBC)
        ! (X, Y): 節点座標, B: ひずみ-変位mat, D: 応力-ひずみmat, TK: 全体剛性mat, E: ヤング率, VNU: ポアソン比, THICK: 厚さ
    double precision :: U(2*MAXNODE), STRAIN(3, 4, MAXELEM), STRESS(3, 4, MAXELEM), XE_C(2)
        ! U: 節点変位ベクトル, STRAIN: 要素のひずみ(epsilon_xx, epsilon_yy, gamma_xy), STRESS: 要素の応力(sigma_xx, sigma_yy, tau_xy)
    double precision :: R1(4), R2(4), DJ(4), DNDR(2, 4, 4), w1(4), w2(4), XE_IP_2(2), XE_IP_3(2)

    open(2, file = 'element_coordinates.txt') ! 節点座標が読み込めてるか確認するための書き出し先ファイル
    open(3, file = '../inputFile/inp_concentrated-f_quadrilateral.txt') ! 節点座標が記述してあるファイル
    open(4, file = 'bmat.txt') ! bmatを記述するファイル
    open(5, file = 'dmat.txt') ! dmatを記述するファイル
    open(6, file = 'kmat.txt') ! 要素剛性マトリックスを記述するファイル
    open(7, file = 'Kmat_diag_P.txt') ! 全体剛性マトリックスの対角項と，B.C.を考慮した場合の右辺ベクトルPを記述するファイル
    open(8, file = 'rst.txt') ! 節点変位ベクトルと各要素の応力
    open(9, file = 'elem_yposi_ystress.csv') ! 中央の要素のy座標（重心位置），y方向垂直応力

    call inp(NNODE, NDOF, NELEM, X, Y, THICK, LNODS, E, VNU, DIM, P, N_BC_GIVEN, I_BC_GIVEN, V_BC_GIVEN) ! 座標系，パラメータのinput
    close(3)

    call mapping(R1, R2, DNDR, W1, W2)
    call dmat(D, E, VNU) ! dmat 作成

    TK(1:NDOF, 1:NDOF) = 0.D0 ! 全体剛性マトリックスを初期化

    do lelem = 1, nelem
        call bmat(lelem, B, DJ, X, Y, LNODS, DNDR) ! bmat 作成(for 4-node)
        call mergekmat(TK, lelem, LNODS, B, DJ, D, THICK, MAXNODE, W1, W2)
    end do

    call bound(TK, P, NDOF, N_BC_GIVEN, I_BC_GIVEN, V_BC_GIVEN, MAXNODE)
    ! do i = 1, NDOF
    !     write(7, '(10F12.4)') TK(i, 1:NDOF)
    ! end do

    write(7, *) 'TK(i, i), P(i)'
    ! write(7, *) TK
    do i = 1, NDOF
        write(7, '(2F12.4)') TK(i, i), P(i)
    end do

    call gauss2(TK, U, P, NDOF, MAXNODE) ! 変位を求める
    ! write(8, *) TK
    write(8, *) 'U(i)'
    ! do inode = 1, NNODE
    !     write(8, '(I5, 2E15.7)') inode, U(2*inode-1), U(2*inode) ! 節点変位を記述
    ! end do
    do i = 1, NDOF
        write(8, '(I5, E15.7)') (i-1)/2 + 1, U(i)
    end do

    do lelem = 1, NELEM
        call bmat(lelem, B, DJ, X, Y, LNODS, DNDR)
        call calc_stress(lelem, STRAIN, STRESS, B, D, U, LNODS) ! ひずみと応力を計算
    end do

    write(8, *) 'STRESS'
    do lelem = 1, NELEM
        write(8, '(I5, 3E15.7)') lelem, 'integ1', STRESS(1:3, 1, lelem) ! 各要素の応力を記述
        write(8, '(I5, 3E15.7)') lelem, 'integ2', STRESS(1:3, 2, lelem)
        write(8, '(I5, 3E15.7)') lelem, 'integ3', STRESS(1:3, 3, lelem)
        write(8, '(I5, 3E15.7)') lelem, 'integ4', STRESS(1:3, 4, lelem)
    end do

    do lelem = 41, 50 ! 中央の要素をピックアップ
        call xe_c_posi(lelem, LNODS, X, Y, U, XE_IP_2, XE_IP_3) ! 中央の要素のy座標（重心位置）をXE_Cに格納
        write(9, *) XE_IP_2(2), ',', STRESS(2, 2, lelem) ! 積分点2のy座標とy方向垂直応力を記述
        write(9, *) XE_IP_3(2), ',', STRESS(2, 3, lelem) ! 積分点3のy座標とy方向垂直応力を記述
    end do
    

    close(2)
    close(4)
    close(5)
    close(6)
    close(7)
    close(8)
    close(9)

    stop
end

! ---------------------------------------------------------------------
!   analytical model input
! ---------------------------------------------------------------------
subroutine inp(nnode, ndof, nelem, x, y, thick, lnods, e, vnu, dim, p, n_bc_given, i_bc_given, v_bc_given)
    implicit none

    integer :: nnode, ndof, nelem, lnods(4, *), dim, n_bc_given, i_bc_given(*), i_load_given, num, idof
    double precision :: thick, x(*), y(*), e, vnu, p(*), v_bc_given(*)

    integer :: inode, lelem
        ! for do文


    ! 機械的特性 ----------
    e = 200D3 ! ヤング率, [MPa]=[N/mm^2], `D3`は10^3を表している
    vnu = 0.4D0 ! ポアソン比, [-]

    thick = 1.D0

    ! inputファイルから読み込み ----------
    read(3, *) nnode, nelem ! 節点数，要素数の読み込み

    do inode = 1, nnode
        read(3, *) x(inode), y(inode) ! 各節点の座標の読み込み
        write(2, *) x(inode), y(inode) ! rst.txtに記述
    end do

    ndof = dim * nnode ! 総自由度数（未知変位成分の数）


    do lelem = 1, nelem
        read(3, *) lnods(1:4, lelem) ! 要素を構成する節点のリスト（必ず反時計回り）
    end do

    ! 境界条件の設定 ----------
    ! 変位境界条件 (x方向)
    num = 0 ! 境界条件の番号，インクリメントする
    idof = 1 ! x方向
    do inode = 1, nnode
        if((x(inode) .eq. 0D0) .and. (y(inode) .eq. 0D0)) then
            num = num + 1
            i_bc_given(num) = inode ! 節点番号
            v_bc_given(num) = 0.D0 ! 強制変位量
            i_bc_given(num) = (i_bc_given(num) - 1) * 2 + idof
        end if
    end do

    ! 変位境界条件（y方向)
    idof = 2 ! y方向
    do inode = 1, nnode
        if(y(inode) .eq. 0D0) then
            num = num + 1
            i_bc_given(num) = inode ! 節点番号
            v_bc_given(num) = 0.D0 ! 強制変位量
            i_bc_given(num) = (i_bc_given(num) - 1) * 2 + idof
        end if 
    end do

    n_bc_given = num

    ! 力学的境界条件
    p(1:ndof) = 0.D0 ! 初期化（力の作用しない自由表面では，節点力はゼロ）
    idof = 2 ! y方向
    do inode = 1, nnode
        if((x(inode) .eq. 0D0) .and. (y(inode) .eq. 50D0)) then
            i_load_given = (inode - 1) * 2 + idof ! 自由度番号
            p(i_load_given) = -10D0 * thick ! 集中荷重: 10N
        end if
    end do


    return
end

! ---------------------------------------------------------------------
!   create B matrix (for each element)
!       - B matrix：shape=(ひずみ成分の数, 要素を構成する節点数x自由度)
! ---------------------------------------------------------------------
subroutine bmat(lelem, b, dj, x, y, lnods, dndr)
    implicit none

    integer :: lelem, lnods(4, *)
        ! lelem: 対象要素の番号, lnods: 対象要素の節点のリスト
    double precision :: x(*), y(*), b(3, 8, 4), dj(4), dndr(2, 4, 4)
        ! (x, y): 各節点の座標のリスト, b: bmat

    integer :: inode, matrix, m, integ, i, j, k, n, i1, i2
    double precision :: xe(2, 4), Ax2, dxdr(2, 2), dxdr_i(2, 2), dndx(2, 4)
        ! xe: 対象要素の節点座標を格納, shape=(2次元平面，要素の頂点の数), Ax2: 要素の面積の2倍

    
    ! 配列xeに対象要素の座標を格納
    do m = 1, 4
        inode = lnods(m, lelem)
        xe(1, m) = x(inode) ! 対象要素のx座標を格納
        xe(2, m) = y(inode) ! 対象要素のy座標を格納
    end do

    do integ = 1, 4
        dxdr(1:2, 1:2) = 0D0
        do i = 1, 2 ! ヤコビマトリクスの作成
            do j = 1, 2
                do k = 1, 4
                dxdr(j, i) = dxdr(j, i) + dndr(j, k, integ) * xe(i, k)
                end do
            end do
        end do

        ! dj(integ) = det[dxdr]を計算する
        dj(integ) = dxdr(1, 1) * dxdr(2, 2) - dxdr(1, 2) * dxdr(2, 1)
        ! dxdrの逆行列をdxdrに記憶する → dr/dx_i になる
        dxdr_i(1, 1) = dxdr(2, 2) / dj(integ)
        dxdr_i(1, 2) = - dxdr(1, 2) / dj(integ)
        dxdr_i(2, 1) = - dxdr(2, 1) / dj(integ)
        dxdr_i(2, 2) = dxdr(1, 1) / dj(integ)

        dndx(1:2, 1:4) = 0.d0 ! 初期化
        do n = 1, 4 
            do i = 1, 2
                do k = 1, 2
                    dndx(i, n) = dndx(i, n) + dxdr_i(i, k) * dndr(k, n, integ) 
                end do
            end do
        end do

        do n = 1, 4 ! 式12.6
            i1 = (n - 1) * 2 + 1
            i2 = (n - 1) * 2 + 2
            b(1, i1, integ) = dndx(1, n)
            b(2, i2, integ) = dndx(2, n)
            b(3, i1, integ) = dndx(2, n)
            b(3, i2, integ) = dndx(1, n)
        end do

    end do ! integ


    do i = 1, 4
        write(4, *) '積分点', i
        do m = 1, 3
            write(4, 100) b(m, 1:8, i) ! bmatをbmat.txtに記述
        end do
    end do
    100 format(8(E12.4, ','))

    return
end

! ---------------------------------------------------------------------
!   create D matrix 
!       - D matrix：shape=(3, 3), 応力-ひずみ間形式の行列
! ---------------------------------------------------------------------
subroutine dmat(d, e, vnu)
    implicit none

    double precision :: d(3, 3), e, vnu
    
    integer :: m
        ! for do文
    double precision :: dCoef


    ! 資料04の式(12)をそのまま
    ! dCoef = e / (1 - vnu**2)
    ! d(1:3, 1:3) = 0.D0 ! 初期化
    ! d(1, 1) = 1
    ! d(1, 2) = vnu
    ! d(2, 1) = vnu
    ! d(2, 2) = 1
    ! d(3, 3) = (1 - vnu) / 2
    ! d = dCoef * d

    ! 資料04の式(6)をべたうち
    dCoef = e / (1 - 2 * vnu) / (1 + vnu)
    d(1:3, 1:3) = 0.D0 ! 初期化
    d(1, 1) = 1 - vnu
    d(1, 2) = vnu
    d(2, 1) = vnu
    d(2, 2) = 1 - vnu
    d(3, 3) = (1 - 2 * vnu) / 2
    d = dCoef * d

    do m = 1, 3
        write(5, 100) d(m, 1:3) ! dmatをdmat.txtに記述
    end do 
    100 format(6(E12.4, ','))

    return
end

! ---------------------------------------------------------------------
!   create K matrix (for each element)
!       - D matrix：shape=(3, 3), 応力-ひずみ間形式の行列
! ---------------------------------------------------------------------
subroutine mergekmat(tk, lelem, lnods, b, dj, d, thick, maxnode, w1, w2)
    implicit none

    integer :: lelem, lnods(4, *), node_no, ip(8), maxnode
        ! lelem: 対象要素の要素番号, lnods: 対象要素の節点のリスト, node_no: 対象節点の番号, 
        ! ip(6): 要素剛性matの行・列と，全体剛性matに組み込む場所の対応関係, db(3, 8): [D][B]を一時記憶
    double precision :: thick, Ax2, Ae, x(maxnode), y(maxnode), xe(2, 4), b(3, 8, 4), d(3, 3), db(3, 8, 4), dk(8, 8, 4), tk(2*maxnode, 2*maxnode)
        ! thick: 厚さ, Ae: 要素の面積, Ax2: Aeの2倍
    double precision :: dj(4), w1(4), w2(4), dk_gi(8, 8)

    integer :: inode, idof, i, j, k, m, integ
        ! for do文


    ! 要素剛性[dk]と全体剛性[tk]の対応関係
    ip(1:8) = 0.D0 ! 初期化
    do inode = 1, 4 ! 要素の4個の節点
        node_no = lnods(inode, lelem) ! 対象節点の番号
        do idof = 1, 2 ! 自由度（degree of freedom）
            ip((inode-1)*2 + idof) = (node_no-1)*2 + idof ! 要素の節点変位ベクトルでの順番 = 全体剛性での自由度番号
        end do
    end do

    ! 配列xeに対象要素の座標を格納
    ! do m = 1, 4
    !     inode = lnods(m, lelem)
    !     xe(1, m) = x(inode) ! 対象要素のx座標を格納
    !     xe(2, m) = y(inode) ! 対象要素のy座標を格納
    ! end do

    dk_gi(1:8, 1:8) = 0.D0
    do integ = 1, 4 ! 積分点ごとに要素剛性行列を計算 (式12.16)
        ! まず，[D(3, 3)][B(3, 8, integ)]を計算し，[db(3, 8, integ)]に記憶
        do i = 1, 3 ! dbの行
            do j = 1, 8 ! db の列
                db(i, j, integ) = 0.D0 ! 初期化
                do k = 1, 3 ! the product of d's row and b's column
                    db(i, j, integ) = db(i, j, integ) + d(i, k) * b(k, j, integ)
                end do
            end do
        end do

        ! 次に，要素剛性mat, [dk(8, 8, integ)] = [B(3, 8, integ)]^T [db(3, 8, integ)]
        do i = 1, 8 ! dk's row
            do j = 1, 8 ! dk's column
                dk(i, j, integ) = 0.D0 ! 初期化
                do k = 1, 3 ! the product of b^T's row and db's column
                    dk(i, j, integ) = dk(i, j, integ) + b(k, i, integ) * db(k, j, integ) ! 転置に注意
                end do
                dk(i, j, integ) = w1(integ) * w2(integ) * dk(i, j, integ) * thick * dj(integ)
                ! ガウス積分を適用
                dk_gi(i, j) = dk_gi(i, j) + dk(i, j, integ) 
            end do
        end do

    end do ! integ

    ! 要素剛性マトリックスを全体剛性マトリックスに組み込む（ミスりがちらしいので注意）
    do i = 1, 8 ! 要素剛性matの行
        do j = 1, 8 ! 要素剛性matの列
            tk(ip(i), ip(j)) = tk(ip(i), ip(j)) + dk_gi(i, j)
        end do 
    end do 

    

    ! 要素番号1--5の要素剛性matを記述したい
    if (lelem .le. 5) then
        write(6, *) '要素番号', lelem, 'の要素剛性'
        write(6, 100) dk_gi(1:8, 1:8)
        100 format(8(E12.4, ','))
    end if

    return
end

! ---------------------------------------------------------------------
!   consider Boundary Condition 
!       - 
! ---------------------------------------------------------------------
subroutine bound(tk, p, ndof, n_bc_given, i_bc_given, v_bc_given, maxnode)
    implicit none

    integer :: ndof, n_bc_given, i_bc_given(*), maxnode
        ! ndof: 総自由度数, n_bc_given: 変位境界条件の数, i_bc_given: 変位境界条件を設定する自由度, maxnode: 最大メモリ？
    double precision :: tk(2*maxnode, 2*maxnode), p(*), v_bc_given(*)
        ! tk: 全体剛性マトリックス, p: 境界条件考慮した時の右辺ベクトル, v_bc_given: 変位境界条件の既知変位の値

    integer :: i, j, K
        ! for do文，jは境界条件を設定する自由度番号を代入する

    ! 元の右辺ベクトルp(*)から系数行列のj列のbarU_j倍をっし引く
    do i = 1, n_bc_given ! i: 境界条件のインデックス
        j = i_bc_given(i) ! 境界条件(i)を設定する自由度番号
        do k = 1, ndof
            p(k) = p(k) - v_bc_given(j) * tk(k, j)
        end do
    end do

    ! 右辺ベクトルp(*)のj行に既知変位の値barU_jを代入する．（上書き）
    do i = 1, n_bc_given
        j = i_bc_given(i) ! 境界条件(i)を設定する自由度番号
        p(j) = v_bc_given(j)
    end do

    ! 系数行列のj行とj列をすべて0にする．対角項(i, j)を1にする
    do i = 1, n_bc_given
        j = i_bc_given(i) ! 境界条件(i)を設定する自由度番号
        do k = 1, ndof
            tk(k, j) = 0
            tk(j, k) = 0
        end do
        tk(j,j) = 1
    end do

end

! ---------------------------------------------------------------------
!   LU decomposition and Gaussian elimination
! ---------------------------------------------------------------------
subroutine gauss2(a, x, b, N, maxnode)

    implicit none
    integer :: N, i, j, k, maxnode
    double precision :: a(2*maxnode, 2*maxnode), x(2*maxnode), b(2*maxnode), sum_au, sum_al, sum_ad, y(N), temp, z(N)

    
    ! 下三角行列 al, 対角行列 ad, 上三角行列 au 成分を左上から右下に向かって逐次的に求めていく ---------------
    j = 1
    ! ad(j, j) = a(j, j)

    j = 2
    a(j-1, j) = a(j-1, j) / a(j-1, j-1)
    a(j, j-1) = a(j, j-1) / a(j-1, j-1)
    a(j, j) = a(j, j) - a(j, j-1) * a(j-1, j-1) * a(j-1, j)

    do j = 3, N ! n: 行列[A]の次元
        a(1, j) = a(1, j) / a(1, 1)
        a(j, 1) = a(j, 1) / a(1, 1)

        do i = 2, j-1
            sum_au = 0.D0 ! 初期化
            do k = 1, i-1
                sum_au = sum_au + a(i, k) * a(k, k) * a(k, j)
                ! write(*, *) a(k,k)
            end do
            
            sum_al = 0.D0 ! 初期化
            do k = 1, i-1
                sum_al = sum_al + a(j, k) * a(k, k) * a(k, i)
            end do
            
            a(i, j) = (a(i, j) - sum_au) / a(i, i)
            a(j, i) = (a(j, i) - sum_al) / a(i, i)
        end do

        sum_ad = 0.D0 ! 初期化
        do k = 1, j-1
            sum_ad = sum_ad + a(j, k) * a(k, k) * a(k, j)
        end do
        a(j, j) = a(j, j) - sum_ad
    end do

    ! 未知ベクトルを求めていく ---------------
    y(1:N) = 0D0 ! 初期化
    ! 前進代入
    y(1) = b(1) ! 1行目
    do i = 2, N
        temp = 0.D0 ! 初期化
        do j = 1, i-1
            temp = temp + a(i, j) * y(j)
        end do  
        y(i) = b(i) - temp
    end do
    
    ! 対角項
    do i = 1, N
        z(i) = y(i) / a(i, i)
    end do

    ! 後退消去
    x(N) = z(N) ! n行目
    do j = N, 2, -1 ! 逆順のインクリメント．3つ目がstep size
        do k = 1, j-1
            z(k) = z(k) - a(k, j) * x(j)
        end do
        x(j-1) = z(j-1)
    end do


    return
end

! ---------------------------------------------------------------------
!   calculate strain and stress
! ---------------------------------------------------------------------
subroutine calc_stress(lelem, strain, stress, b, d, u, lnods)
    implicit none

    integer :: lelem, lnods(4, *), m, inode, integ
    double precision :: strain(3, 4, *), stress(3, 4, *), b(3, 8, 4), d(3, 3), u(*), ue(8)
    ! 対象要素の節点変位を持ってくる
    ue(1:8) = 0.D0 ! 初期化
    do m = 1, 4
        inode = lnods(m, lelem)
        ue(2*m-1) = u(2*inode-1)
        ue(2*m) = u(2*inode)
    end do

    do integ = 1, 4
        ! ひずみ-変位関係式でひずみを計算
        strain(1:3, integ, lelem) = matmul(b(1:3, 1:8, integ), ue)
        ! 応力-ひずみ間形式で応力を計算
        stress(1:3, integ, lelem) = matmul(d, matmul(b(1:3, 1:8, integ), ue))
    end do

    return 
end

! ---------------------------------------------------------------------
!   変形後の重心の座標を計算
! ---------------------------------------------------------------------
subroutine xe_c_posi(lelem, lnods, x, y, u, xe_c)
    implicit none

    integer :: lelem, lnods(4, *), m, inode
    double precision :: x(*), y(*), u(*), xe(2, 3), xe_c(2)
    ! 配列xeに対象要素の変形後の座標を格納
    do m = 1, 4
        inode = lnods(m, lelem)
        xe(1, m) = x(inode) + u(2*inode-1)
        xe(2, m) = y(inode) + u(2*inode)
    end do
    ! 重心の座標
    xe_c(1) = (xe(1, 1) + xe(1, 2) + xe(1, 3)) / 3
    xe_c(2) = (xe(2, 1) + xe(2, 2) + xe(2, 3)) / 3

    return 
end

! ---------------------------------------------------------------------
!   isoparametric mapping 
! ---------------------------------------------------------------------
subroutine mapping(r1, r2, dndr, w1, w2)
    implicit none

    integer :: integ
    double precision :: r1(4), r2(4), dndr(2, 4, 4), w1(4), w2(4)

    r1(1) = -1.D0
    r1(2) = 1.D0
    r1(3) = 1.D0
    r1(4) = -1.D0
    r1(1:4) = r1(1:4) / dsqrt(3D0)
    w1(1:4) = 1D0

    r2(1) = -1.D0
    r2(2) = -1.D0
    r2(3) = 1.D0
    r2(4) = 1.D0
    r2(1:4) = r2(1:4) / dsqrt(3D0)
    w2(1:4) = 1D0

    do integ = 1, 4
        dndr(1, 1, integ) = -0.25d0 * (1.d0 - r2(integ))
        dndr(2, 1, integ) = -0.25d0 * (1.d0 - r1(integ))
        dndr(1, 2, integ) = 0.25d0 * (1.d0 - r2(integ))
        dndr(2, 2, integ) = -0.25d0 * (1.d0 + r1(integ))
        dndr(1, 3, integ) = 0.25d0 * (1.d0 + r2(integ))
        dndr(2, 3, integ) = 0.25d0 * (1.d0 + r1(integ))
        dndr(1, 4, integ) = -0.25d0 * (1.d0 + r2(integ))
        dndr(2, 4, integ) = 0.25d0 * (1.d0 - r1(integ))
    end do
    
    return
end

! ---------------------------------------------------------------------
!   積分点（2, 3）の座標を計算
! ---------------------------------------------------------------------
subroutine xe_ip_posi(lelem, lnods, x, y, u, xe_ip_2, xe_ip_3)
    implicit none

    integer :: lelem, lnods(4, *), m, inode
    double precision :: x(*), y(*), u(*), xe(2, 3), xe_ip_2(2), xe_ip_3(2)
    ! 配列xeに対象要素の変形後の座標を格納
    do m = 1, 4
        inode = lnods(m, lelem)
        xe(1, m) = x(inode) + u(2*inode-1)
        xe(2, m) = y(inode) + u(2*inode)
    end do
    ! 積分点の座標
    xe_ip_2(1) = (xe(1, 1) + xe(1, 2) + xe(1, 3)) / 3 + dsqrt(3D0)
    xe_ip_2(2) = (xe(2, 1) + xe(2, 2) + xe(2, 3)) / 3 - dsqrt(3D0)
    
    xe_ip_3(1) = (xe(1, 1) + xe(1, 2) + xe(1, 3)) / 3 + dsqrt(3D0)
    xe_ip_3(2) = (xe(2, 1) + xe(2, 2) + xe(2, 3)) / 3 + dsqrt(3D0)

    

    return 
end




