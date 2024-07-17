!
!  Copyright (C) 2018 by the authors of the RAYLEIGH code.
!
!  This file is part of RAYLEIGH.
!
!  RAYLEIGH is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3, or (at your option)
!  any later version.
!
!  RAYLEIGH is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with RAYLEIGH; see the file LICENSE.  If not see
!  <http://www.gnu.org/licenses/>.
!

Module Sphere_Linear_Terms
    Use Load_Balance, Only : mp_lm_values, l_lm_values, my_num_lm, m_lm_values, my_lm_min, my_nl_lm, my_nm_lm, my_lm_lval, my_lm_max
    Use Controls
    Use ProblemSize
    Use Linear_Solve
    Use Fields
    Use BoundaryConditions
    Use Timers
    Use ClockInfo
    Use PDE_Coefficients
    Use Math_Constants
    Implicit None
    Real*8, Allocatable :: Lconservation_weights(:)

Contains
    Subroutine Linear_Init()
        Implicit None
        Real*8 :: amp, T,arg
        Integer :: n, r, m, nm
        !Depending on process layout, some ranks may not participate in the solve
        If (my_num_lm .gt. 0) Then 
 
            Call Initialize_Linear_System()

            If (strict_L_conservation) Then

                Allocate(Lconservation_weights(1:N_R))
                Lconservation_weights(1:N_R) = 0.0d0
                nm = 0
                Do m = 1, gridcp%domain_count
                    Do n = 1, gridcp%npoly(m)
                        Do r = 1, gridcp%npoly(m)
                            T = gridcp%dcheby(m)%data(r,n,0)
                            Lconservation_weights(n+nm) = Lconservation_weights(n+nm) + radial_integral_weights(r+nm) * T
                        Enddo
                    Enddo
                    Lconservation_weights( nm+(2*gridcp%npoly(m))/3+1:nm+gridcp%npoly(m) ) = 0.0d0  ! De-Alias
                    nm = nm + gridcp%npoly(m)
                Enddo
            Endif
        Endif
    End Subroutine Linear_Init

    Subroutine Reset_Linear_Equations()
        Implicit None
        Real*8 :: rhs_factor, lhs_factor
        Call Reset_Equation_Coefficients()    ! Zero out all coefficients
        lhs_factor = -deltat*alpha_implicit    ! Crank Nicolson scheme - alpha = 0.5
        rhs_factor = deltat*(1.0d0-alpha_implicit)
        Call Set_Time_Factors(lhs_factor,rhs_factor)
        Call Load_Linear_Coefficients()
        Call LU_Decompose_Matrices()
    End Subroutine Reset_Linear_Equations


    Subroutine Initialize_Linear_System
        Implicit None
        Integer :: lp, l, nlinks
        Integer, Allocatable :: eq_links(:), var_links(:)
        Type(Cheby_Grid), Pointer :: gridpointer
        nullify(gridpointer)
        gridpointer => gridcp
        If (chebyshev) Call Use_Chebyshev(gridpointer)    ! Turns chebyshev mode to "on" for the linear solve

        Call Initialize_Equation_Set(n_equations,n_variables,N_R,my_nl_lm, my_nm_lm,2)

        Do lp = 1, my_nl_lm
            l = my_lm_lval(lp)
            Call Initialize_Equation_Coefficients(vreq ,     vr, 1, lp)
            Call Initialize_Equation_Coefficients(vteq , vtheta, 1, lp)
            Call Initialize_Equation_Coefficients(vpeq ,   vphi, 1, lp)
            Call Initialize_Equation_Coefficients(teq  ,   tvar, 1, lp)
            Call Initialize_Equation_Coefficients(rhoeq, rhovar, 1, lp)

            If (l .ne. 0) Then
                If (magnetism) Then
                    Call Initialize_Equation_Coefficients(ceq,cvar, 2,lp)
                    Call Initialize_Equation_Coefficients(aeq,avar, 2,lp)
                Endif
            Endif
        Enddo
        Call Finalize_Equations()
        If (bandsolve) Call Use_BandSolve()
        If (sparsesolve) Call Use_SparseSolve()

    End Subroutine Initialize_Linear_System

    Subroutine Load_Linear_Coefficients()
        Implicit None

        Real*8, Allocatable :: H_Laplacian(:), amp(:)
        Integer :: l, lp
        Real*8 :: diff_factor,ell_term
        !rmin_norm
        diff_factor = 1.0d0 ! hyperdiffusion factor (if desired, 1.0d0 is equivalent to no hyperdiffusion)
        Allocate(amp(1:N_R))
        Allocate(H_Laplacian(1:N_R))
        Do lp = 1, my_nl_lm
            If (bandsolve) Call DeAllocate_LHS(lp)
            Call Allocate_LHS(lp)
            l = my_lm_lval(lp)

            If (hyperdiffusion) Then
                ell_term = ((l-1.0d0)/(l_max-1.0d0))**hyperdiffusion_beta
                diff_factor = 1.0d0+hyperdiffusion_alpha*ell_term
            Endif
            H_Laplacian = - l_l_plus1(l) * OneOverRSquared

            !==================================================
            
            amp = 1.0d0
            Call add_implicit_term( vreq,     vr, 0, amp,lp,static = .true.)
            Call add_implicit_term( vteq, vtheta, 0, amp,lp,static = .true.)
            Call add_implicit_term( vpeq,   vphi, 0, amp,lp,static = .true.)
            Call add_implicit_term(  teq,   tvar, 0, amp,lp,static = .true.)
            Call add_implicit_term(rhoeq, rhovar, 0, amp,lp,static = .true.)


            If (l .ne. 0) Then
               
                If (magnetism) Then
                    !=========================================
                    !  Btor Equation
                    amp = 1.0d0
                    Call add_implicit_term(aeq,avar, 0, amp,lp, static = .true.)    ! Time-independent piece

                    amp = H_Laplacian*eta*diff_factor
                    Call add_implicit_term(aeq,avar, 0, amp,lp)

                    amp = 1.0d0*eta*diff_factor
                    Call add_implicit_term(aeq,avar, 2, amp,lp)

                    ! Eta variation in radius
                    amp = A_Diffusion_Coefs_1*diff_factor
                    Call add_implicit_term(aeq,avar,1,amp,lp)

                    !=========================================
                    !  Bpol Equation
                    amp = 1.0d0
                    Call add_implicit_term(ceq,cvar, 0, amp,lp, static = .true.)    ! Time-independent piece

                    amp = H_Laplacian*eta*diff_factor
                    Call add_implicit_term(ceq,cvar, 0, amp,lp)

                    amp = 1.0d0*eta*diff_factor
                    Call add_implicit_term(ceq,cvar, 2, amp,lp)
                Endif

                ! If band solve, do the redefinition of the matrix here

            Endif
            Call Set_Boundary_Conditions(lp)
            If (sparsesolve) Then
                !Write(6,*)'matrix: ', weq,lp, my_rank, l
                Call Sparse_Load(weq,lp)
                !Write(6,*)'matrix: ', zeq,lp,my_rank, l
                Call Sparse_Load(zeq,lp)
                If (magnetism) Then
                    Call Sparse_Load(aeq,lp)
                    Call Sparse_Load(ceq,lp)
                Endif
            Endif


            If (bandsolve) Then
                Call Band_Arrange(weq,lp)
                Call Band_Arrange(zeq,lp)
                If (magnetism) Then
                    Call Band_Arrange(aeq,lp)
                    Call Band_Arrange(ceq,lp)
                Endif
            Endif

        Enddo
        DeAllocate(amp)
        DeAllocate(H_Laplacian)
    End Subroutine Load_Linear_Coefficients

    Subroutine Set_Boundary_Conditions(mode_ind)
        ! Modified version of set_boundary_conditions
        ! Designed to work with more memory friendly logic
        ! Sets boundary condition of indicated l-value (index lp)
        ! only.  Does not loop over lp.
        Implicit None
        Real*8 :: samp,one
        Integer, Intent(In) :: mode_ind
        Integer :: l, r,lp
        one = 1.0d0
        lp = mode_ind

        l = my_lm_lval(lp)

        !*******************************************************
        !        Clear the boundary rows
        Call Clear_Row(vreq,lp,1)
        Call Clear_Row(vreq,lp,N_R)
        Call Clear_Row(vteq,lp,1)
        Call Clear_Row(vteq,lp,N_R)
        Call Clear_Row(vpeq,lp,1)
        Call Clear_Row(vpeq,lp,N_R)
        Call Clear_Row(teq,lp,1)
        Call Clear_Row(teq,lp,N_R)

        ! "1" denotes linking at index 1, starting in domain 2
        ! "2" denotes linking at index npoly, starting in domain 1


        ! For each variable except for rho, var and var' are continuous

        Call FEContinuity(vreq,lp,vr,2,0)     ! var
        Call FEContinuity(vreq,lp,vr,1,1)     ! var' 

        Call FEContinuity(vteq,lp,vtheta,2,0) ! var
        Call FEContinuity(vteq,lp,vtheta,1,1) ! var' 

        Call FEContinuity(vpeq,lp,vphi,2,0)   ! var
        Call FEContinuity(vpeq,lp,vphi,1,1)   ! var' 


        !Call FEContinuity(teq,lp,tvar,2,0)   ! var
        !Call FEContinuity(teq,lp,tvar,1,1)   ! var' 

        Call FEContinuity(rhoeq,lp,rhovar,1,0)   ! rho is continous

        !***********************************************************
        ! Temperature Boundary Conditions 
        r = 1
        If (fix_tvar_top) Then

            Call Load_BC(lp,r,teq,tvar,one,0)    !upper boundary

        Endif
        If (fix_dtdr_top) Then
            Call Load_BC(lp,r,teq,tvar,one,1)
        Endif

        r = N_R
        If (fix_tvar_bottom) Then
            Call Load_BC(lp,r,teq,tvar,one,0)    ! lower boundary
        Endif
        If (fix_dtdr_bottom) Then
            Call Load_BC(lp,r,teq,tvar,one,1)
        Endif

        !************************************************************
        ! Velocity Boundary Conditions

        ! Impenetrable top and bottom
        ! vr vanishes at the boundaries
        r = 1
        Call Load_BC(lp,r,vreq,vr,one,0)
        r = N_R
        Call Load_BC(lp,r,vreq,vr,one,0)

        ! Density
        ! I'm a little uncertain what a meaningful boundary condition is here
        ! Let's try density gradient = 0 at upper boundary
        ! Possibly enforce mass conservation through integral constraint on
        ! ell = 0?

        !r = 1
        !Call Load_BC(lp,r,rhoeq,rhovar,one,1)

        If (no_slip_top) Then
            r = 1

            Call Load_BC(lp,r,vteq,vtheta,one,0)
            Call Load_BC(lp,r,vpeq,vphi,one,0)
        Else
            ! Else stress-free
            r = 1
            samp = -1.0d0/radius(r)
            Call Load_BC(lp,r,vteq,vtheta,one,1)
            Call Load_BC(lp,r,vteq,vtheta,samp,0)


            Call Load_BC(lp,r,vpeq,vphi,one,1)
            Call Load_BC(lp,r,vpeq,vphi,samp,0)
        Endif

        If (no_slip_bottom) Then
            r = N_R
            Call Load_BC(lp,r,vteq,vtheta,one,0)
            Call Load_BC(lp,r,vpeq,vphi,one,0)

        Else
            !stress_free_bottom
            r = N_R
            samp = -1.0d0/radius(r)
            Call Load_BC(lp,r,vteq,vtheta,one,1)
            Call Load_BC(lp,r,vteq,vtheta,samp,0)
            Call Load_BC(lp,r,vpeq,vphi,one,1)
            Call Load_BC(lp,r,vpeq,vphi,samp,0)
        Endif



        If (l .ne. 0) Then

            !*******************************************************
            !        Magnetic Boundary Conditions

            If (Magnetism) Then
                !  Clear the boundary rows
                Call Clear_Row(ceq,lp,1)
                Call Clear_Row(ceq,lp,N_R)
                Call Clear_Row(aeq,lp,1)
                Call Clear_Row(aeq,lp,N_R)


                Call FEContinuity(aeq,lp,avar,2,0)   ! A is continuous
                Call FEContinuity(aeq,lp,avar,1,1)          ! A' is continuous

                Call FEContinuity(ceq,lp,cvar,2,0)   ! C is continuous
                Call FEContinuity(ceq,lp,cvar,1,1)          ! C' is continuous


                ! Match to a potential field at top and bottom
                ! Btor = 0 at top and bottom
                r = 1
                Call Load_BC(lp,r,aeq,avar,one,0)
                r = N_R
                Call Load_BC(lp,r,aeq,avar,one,0)

                ! dBpol/dr+ell*Bpol/r = 0 at outer boundary
                r = 1
                Call Load_BC(lp,r,ceq,cvar,one,1)
                samp = my_lm_lval(lp)*one_over_r(r)
                Call Load_BC(lp,r,ceq,cvar,samp,0)

                ! dBpol/dr-(ell+1)*Bpol/r = 0 at inner boundary
                r = N_R
                Call Load_BC(lp,r,ceq,cvar,one,1)
                samp = - (l+1)*One_Over_R(r)
                Call Load_BC(lp,r,ceq,cvar,samp,0)


                If (fix_poloidalfield_top) Then
                    Call Clear_Row(ceq,lp,1)
                    Call Clear_Row(aeq,lp,1)
                    r = 1
                    Call Load_BC(lp,r,ceq,cvar,one,0)
                    Call Load_BC(lp,r,aeq,avar,one,0)
                Endif
                If (fix_poloidalfield_bottom) Then
                    Call Clear_Row(aeq,lp,N_R)
                    Call Clear_Row(ceq,lp,N_R)
                    r = N_R
                    Call Load_BC(lp,r,ceq,cvar,one,0)
                    Call Load_BC(lp,r,aeq,avar,one,0)
                Endif

            Endif    ! Magnetism


        Endif ! l = 0 or not


    End Subroutine Set_Boundary_Conditions

    Subroutine Enforce_Boundary_Conditions
        Implicit None
        Real*8 :: bc_val
        Integer :: uind, lind
        Integer :: real_ind, imag_ind

        Call Apply_Boundary_Mask(bc_values,vreq)
        Call Apply_Boundary_Mask(bc_values,vteq)
        Call Apply_Boundary_Mask(bc_values,vpeq)
        Call Apply_Boundary_Mask(bc_values,teq)
        Call Domain_Continuity()

    End Subroutine Enforce_Boundary_Conditions

    !//////////////////////////////////////////////////////////////////////////////
    ! The domain continuity routine still needs to be updated for compressible mode
    ! Unless we run with multiple chebyshev domains, however, this isn't an issue
    Subroutine Domain_Continuity()
        Implicit None
        Integer :: l, indx, ii,lp, j, k,n
        ! start applying the boundary and continuity conditions by setting
        ! the appropriate right hand sides.

        ! Will wrap this more into Linear_Solve.F90 soon
        ii = 2*N_R


        indx = 1

        Do lp = 1, my_nl_lm
            n_m = my_nm_lm(lp)-1    ! really n_m, but the indexing below is from the old implicit solv
            l = my_lm_lval(lp)

            If (chebyshev) Then
                ! Apply continuity conditions across the subdomains

                j = 0
                do n = 1, ndomains-1
                    j = j+gridcp%npoly(n)
                    equation_set(1,weq)%RHS(      j, : , indx:indx+n_m) = 0.0d0
                    equation_set(1,weq)%RHS(2*N_R+j, : , indx:indx+n_m) = 0.0d0
                    if (l .ne. 0) then
                        equation_set(1,zeq)%RHS(      j, : , indx:indx+n_m) = 0.0d0
                        equation_set(1,weq)%RHS(N_R  +j, : , indx:indx+n_m) = 0.0d0 ! dpdr continuity - only for ell =/ 0
                    endif
                Enddo


                j = 1
                Do n = 1, ndomains-1
                    j = j+gridcp%npoly(n)
                    equation_set(1,weq)%RHS(      j  ,:, indx:indx+n_m) = 0.0d0
                    equation_set(1,weq)%RHS(  N_R+j  ,:, indx:indx+n_m) = 0.0d0
                    equation_set(1,weq)%RHS(2*N_R+j  ,:, indx:indx+n_m) = 0.0d0
                    if (l .ne. 0) equation_set(1,zeq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
                Enddo

                If (Magnetism) Then
                    if (l .ne. 0) then
                        j = 0
                        Do n = 1, ndomains-1
                            j = j+gridcp%npoly(n)
                            equation_set(1,aeq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
                            equation_set(1,ceq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
                        Enddo

                        j = 1
                        Do n = 1, ndomains-1
                            j = j+gridcp%npoly(n)
                            equation_set(1,aeq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
                            equation_set(1,ceq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
                        Enddo
                    endif
                Endif


            Endif

            indx = indx + n_m + 1
        Enddo
    End Subroutine Domain_Continuity

End Module Sphere_Linear_Terms
