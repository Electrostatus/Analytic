module analytic  ! DESCENDING ORDER!
    implicit none
    ! A collection of general purpose analytic formulas
    ! for polynomials of degree 0 through 4 

    ! Polynomials with degrees higher than four do not have analytic solutions

    ! this file accepts the coefficients in descending order: cof(1) /= 0
    !    (cof(1) * x^n + cof(2) * x^(n-1) + cof(3) * x^(n-2) + ...)
    ! for coefficients in ascending order, see 'analytic_asc.f90'

    integer, parameter  :: dp = selected_real_kind(8)  ! default precision
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

    ! second and third cubic unity roots (first is just 1)
    real(dp), parameter :: d60 = sin(pi / 3.0_dp), nd60 = -sin(pi / 3.0_dp)
    complex(dp), parameter :: cu2 = (-0.5_dp, d60), cu3 = (-0.5_dp, nd60)

    private
    public :: root_0, root_1, root_2, root_3, root_4

    contains
        ! cof - complex array of coefficients,
        !       in descending order; cof(1) /= 0
        ! rts - complex array of found roots
        subroutine root_0(cof, rts)  ! constant,  a = 0
            complex(dp), intent(in)  :: cof(:)
            complex(dp), intent(out) :: rts(:)
            rts(1) = 0.0_dp
        end subroutine root_0

        subroutine root_1(cof, rts)  ! linear,    ax + b = 0
            complex(dp), intent(in)  :: cof(:)
            complex(dp), intent(out) :: rts(:)
            rts(1) = -cof(1) / cof(2)
        end subroutine root_1

        subroutine root_2(cof, rts)  ! quadratic, ax^2 + bx + c = 0
            complex(dp), intent(in)  :: cof(:)
            complex(dp), intent(out) :: rts(:)
            complex(dp) :: p, q

            p = sqrt(cof(2) * cof(2) - 4.0_dp * cof(1) * cof(3))
            q = -2.0_dp * cof(1)
            rts(1) = (cof(2) - p) / q
            rts(2) = (cof(2) + p) / q
        end subroutine root_2

        subroutine root_3(cof, rts)  ! cubic,     ax^3 + bx^2 + cx + d = 0
            complex(dp) :: abc, bbb, aad, dd, d0, d1, cin, cc, p, x1
            complex(dp), intent(in)  :: cof(:)
            complex(dp), intent(out) :: rts(:)

            abc = cof(1) * cof(2) * cof(3)
            bbb = cof(2) * cof(2) * cof(2)
            aad = cof(1) * cof(1) * cof(4)

            dd = 18.0_dp * abc * cof(4) - 4.0_dp * bbb * cof(4)
            dd = dd + cof(2) * cof(2) * cof(3) * cof(3) - 4.0_dp * cof(1) * cof(3) * cof(3) * cof(3)
            dd = dd - 27.0_dp * aad * cof(4)
            d0 = cof(2) * cof(2) - 3.0_dp * cof(1) * cof(3)

            if ((dd == 0).and.(d0 == 0)) then  ! triple root
                x1 = -cof(2) / (3.0_dp * cof(1))
                rts(1) = x1
                rts(2) = x1
                rts(3) = x1
            else if ((dd == 0).and.(d0 /= 0)) then  ! double root, simple root
                x1 = ((9.0_dp * cof(1) * cof(4) - cof(2) * cof(3)) / (2.0_dp * d0))
                rts(1) = x1
                rts(2) = x1
                rts(3) = (4.0_dp * abc - 9.0_dp * aad - bbb) / (cof(1) * d0)
            else
                d1 = 2.0_dp * bbb - 9.0_dp * abc + 27.0_dp * aad

                if (d0 == 0) then
                    cin = d1
                else
                    cin = (d1 - sqrt(-27.0_dp * cof(1) * cof(1) * dd)) / 2.0_dp
                end if
                cc = cin ** (1.0_dp / 3.0_dp)
                p = (-1.0_dp / (3.0_dp * cof(1)))

                rts(1) = p * (cof(2) + cc + d0 / cc)
                rts(2) = p * (cof(2) + cu2 * cc + d0 / (cu2 * cc))
                rts(3) = p * (cof(2) + cu3 * cc + d0 / (cu3 * cc))
            end if
        end subroutine root_3

        subroutine root_4(cof, rts)  ! quartic,   ax^4 + bx^3 + cx^2 + dx + e = 0
            complex(dp) :: a, b, c, d, a2, b2, bq, c1
            complex(dp) :: usq, usq2, u4, blkp, blkm
            complex(dp) :: p1, v, u, uu, v3, p2, p3
            complex(dp), intent(in)  :: cof(:)
            complex(dp), intent(out) :: rts(:)

            a = cof(2) / cof(1)
            b = cof(3) / cof(1)
            c = cof(4) / cof(1)
            d = cof(5) / cof(1)

            a2 = a * a; b2 = b * b

            bq = -(2.0_dp * b2 * b) + 9.0_dp * a * b * c
            bq = bq - 27.0_dp * (c * c + a2 * d) + 72.0_dp * b * d
            c1 = b2 - 3.0_dp * a * c + 12.0_dp * d

            p1 = sqrt(bq * bq - (4.0_dp * c1 * c1 * c1))
            v = -(bq - p1) / 2.0_dp  ! solve quadratic
            if (v == 0) v = -(bq + p1) / 2.0_dp  ! choose non zero quad root

            u = a2 / 4.0_dp - (2.0_dp * b) / 3.0_dp
            if (v == 0) then  ! both quad roots zero, uu simplifies to u
                uu = u
            else
                v3 = (v ** (1.0_dp / 3.0_dp)) * cu2
                uu = u + (1.0_dp / 3.0_dp) * (v3 + c1 / v3)
            end if

            p1 = - a / 4.0_dp

            if (uu == 0) then  ! degenerate, quadruple root
                rts(1) = p1
                rts(2) = p1
                rts(3) = p1
                rts(4) = p1
            else
                p2 = 3.0_dp * a2 - 8.0_dp * b - 4.0_dp * uu
                p3 = -(a2 * a) + 4.0_dp * a * b - 8.0_dp * c

                usq = sqrt(uu)
                usq2 = usq / 2.0_dp
                u4 = uu / 4.0_dp

                blkp = (1.0_dp / 4.0_dp) * sqrt(p2 + p3 / usq)
                blkm = (1.0_dp / 4.0_dp) * sqrt(p2 - p3 / usq)

                rts(1) = p1 + usq2 + blkp
                rts(2) = p1 - usq2 + blkm
                rts(3) = p1 + usq2 - blkp
                rts(4) = p1 - usq2 - blkm
            end if
        end subroutine root_4

end module analytic

!program example
!use analytic  ! descending order
!complex(8) :: cofs(5), rts(4)
!
!cofs = [3, 0, 0, 0, 0]
!call root_0(cofs, rts)
!print *,'root_0 in : ',cofs(1:1)
!print *,'root_0 out: ', rts(1:1)
!print *, ''
!
!cofs = [3, 2, 0, 0, 0]
!call root_1(cofs, rts)
!print *,'root_1 in : ',cofs(1:2)
!print *,'root_1 out: ', rts(1:1)
!print *, ''
!
!cofs = [3, 2, 7, 0, 0]
!call root_2(cofs, rts)
!print *,'root_2 in : ',cofs(1:3)
!print *,'root_2 out: ', rts(1:2)
!print *, ''
!
!cofs = [3, 2, 7, 13, 0]
!call root_3(cofs, rts)
!print *,'root_3 in : ',cofs(1:4)
!print *,'root_3 out: ', rts(1:3)
!print *, ''
!
!cofs = [3, 2, 7, 13, 5]
!call root_4(cofs, rts)
!print *,'root_4 in : ',cofs(1:5)
!print *,'root_4 out: ', rts(1:4)
!end program example