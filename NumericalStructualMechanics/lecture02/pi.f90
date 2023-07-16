program pi
! 円周率をTaylor展開から算出するプログラム

    integer :: a = 1e6
    double precision result
    double precision, external :: invtan

    result = 4.0d0 * inv_tan(1.0d0, a)
    
    open(35, file='pi.txt')
    write(35, *) "the result of pi.f90: ", result
    write(35, *) "the number of terms: ", a
    close(35)

    stop
contains
    ! inverse tangentの関数
    double precision function inv_tan(x, a)
        double precision x ! 引数
        integer a ! 引数
        double precision :: ans=0.0d0 
        integer i

        do i = 1, a
            ans = ans + (-1.0d0)**(i + 1) * x**(2*i - 1) / (2*i - 1)
        enddo

        inv_tan = ans

    end function inv_tan

end program 


