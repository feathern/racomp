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
    Use ReferenceState
    Use TransportCoefficients
    Use Math_Constants
    Implicit None
    Real*8, Allocatable :: Lconservation_weights(:)

Contains
    Subroutine Linear_Init()
        Implicit None
        Real*8 :: amp, T,arg
        Integer :: n, r
        !Depending on process layout, some ranks may not participate in the solve
        If (my_num_lm .gt. 0) Then

            Call Initialize_Linear_System()

            If (strict_L_conservation) Then

                Allocate(Lconservation_weights(1:N_R))
                Lconservation_weights(1:N_R) = 0.0d0

                Do n = 1, N_R
                    Do r = 1, N_R
                        T = gridcp%dcheby(1)%data(r,n,0)
                        Lconservation_weights(n) = Lconservation_weights(n) + radial_integral_weights(r) * T
                    Enddo
                Enddo

                Lconservation_weights( (2*N_R)/3+1: ) = 0.0d0  ! De-Alias

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
        Integer :: neq, nvar,lp, l, nlinks
        Integer, Allocatable :: eq_links(:), var_links(:)
        Type(Cheby_Grid), Pointer :: gridpointer
        If (magnetism) Then
            neq  = 6
            nvar = 6
        Else
            neq  = 5
            nvar = 5
        Endif
        nullify(gridpointer)
        gridpointer => gridcp
        If (chebyshev) Call Use_Chebyshev(gridpointer)    ! Turns chebyshev mode to "on" for the linear solve

        Call Initialize_Equation_Set(neq,nvar,N_R,my_nl_lm, my_nm_lm,2)

        Do lp = 1, my_nl_lm
            l = my_lm_lval(lp)

            ! W equation
            Call Initialize_Equation_Coefficients(vreq ,     vr, 1, lp)
            Call Initialize_Equation_Coefficients(vteq , vtheta, 1, lp)
            Call Initialize_Equation_Coefficients(vpeq ,   vphi, 1, lp)
            Call Initialize_Equation_Coefficients(teq  ,   tvar, 1, lp)
            Call Initialize_Equation_Coefficients(rhoeq, rhovar, 1, lp)

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
            Call Clear_Row(rhoeq,lp,1)  ! Density gets one B.C., I think... (NF, Jun 18, 2019)

            ! "1" denotes linking at index 1, starting in domain 2
            ! "2" denotes linking at index npoly, starting in domain 1


            ! For each variable except for rho, var and var' are continuous

            Call FEContinuity(vreq,lp,vr,2,0)     ! var
            Call FEContinuity(vreq,lp,vr,1,1)     ! var' 

            Call FEContinuity(vteq,lp,vtheta,2,0) ! var
            Call FEContinuity(vteq,lp,vtheta,1,1) ! var' 

            Call FEContinuity(vpeq,lp,vphi,2,0)   ! var
            Call FEContinuity(vpeq,lp,vphi,1,1)   ! var' 


            Call FEContinuity(teq,lp,tvar,2,0)   ! var
            Call FEContinuity(teq,lp,tvar,1,1)   ! var' 

            Call FEContinuity(rhoeq,lp,rhovar,1,0)   ! rho is continous




            ! Temperature Boundary Conditions 
            r = 1
            If (fix_tvar_top) Then
                If (.not. fix_divrfc_top) Then
                    Call Load_BC(lp,r,teq,tvar,one,0)    !upper boundary
                Endif
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

            !///////////////////////////////////////////////////////////////////
            ! These four different boundary conditions are similar, though
            ! slightly different, in nature.  The idea is to try some boundary conditions
            ! that allow entropy and it's derivatives to vary on the boundary
            ! Either Del dot Grad S, Del dot F_conductive, or Del_r dot Grad S or F_conductive is zero






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

            r = 1
            Call Load_BC(lp,r,rhoeq,rhovar,one,1)

            If (no_slip_top) Then
                r = 1

                Call Load_BC(lp,r,vteq,vtheta,one,0)
                Call Load_BC(lp,r,vpeq,vphi,one,0)
            Else
                ! Else stress-free
                r = 1
                samp = -2.0d0/radius(r)
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
                samp = -2.0d0/radius(r)
                Call Load_BC(lp,r,vteq,vtheta,one,1)
                Call Load_BC(lp,r,vteq,vtheta,samp,0)
                Call Load_BC(lp,r,vpeq,vphi,one,1)
                Call Load_BC(lp,r,vpeq,vphi,samp,0)
            Endif

            If ((l .eq. 1) .and. (strict_L_Conservation) ) then
                Call Clear_Row(zeq,lp,1)
                Call Load_BC(lp,1,zeq,zvar,one,0,integral = Lconservation_weights)
            Endif





    End Subroutine Set_Boundary_Conditions



    Subroutine Enforce_Boundary_Conditions()
        Implicit None
        Integer :: l, indx, ii,lp, j, k,n
        ! start applying the boundary and continuity conditions by setting
        ! the appropriate right hand sides.

        ! This is ugly, and I have no idea how to make this pretty.
        ! Might zero these by default and then call a routine for exceptions
        ! such as fixed entropy top.
        ii = 2*N_R


        indx = 1

        Do lp = 1, my_nl_lm
            n_m = my_nm_lm(lp)-1    ! really n_m, but the indexing below is from the old implicit solv
            l = my_lm_lval(lp)



            equation_set(1,vreq)%RHS(1  ,:,indx:indx+n_m)   = Zero
            equation_set(1,vreq)%RHS(N_R,:,indx:indx+n_m)   = Zero

            equation_set(1,vteq)%RHS(1  ,:,indx:indx+n_m)   = Zero
            equation_set(1,vteq)%RHS(N_R,:,indx:indx+n_m)   = Zero

            equation_set(1,vpeq)%RHS(1  ,:,indx:indx+n_m)   = Zero
            equation_set(1,vpeq)%RHS(N_R,:,indx:indx+n_m)   = Zero

            equation_set(1,teq)%RHS(1  ,:,indx:indx+n_m)    = Zero
            equation_set(1,teq)%RHS(N_R,:,indx:indx+n_m)    = Zero

            equation_set(1,rhoeq)%RHS(1  ,:,indx:indx+n_m)  = Zero

            If (l .eq. 0) Then
                equation_set(1, vreq)%RHS(:,2,indx)   = Zero ! no imaginary part for any ell=0 equations
                equation_set(1, vteq)%RHS(:,2,indx)   = Zero
                equation_set(1, vpeq)%RHS(:,2,indx)   = Zero
                equation_set(1,  teq)%RHS(:,2,indx)   = Zero
                equation_set(1,rhoeq)%RHS(:,2,indx)   = Zero

                If (fix_tvar_top) Then
                    !Top temperature (in spectral space, but BC's specified in physical space
                    !    so multiply by sqrt4pi)
                    equation_set(1,teq)%RHS(1,1,indx)   = T_Top*sqrt(4.0D0*Pi)
                Endif
                If (fix_dtdr_top) Then
                    equation_set(1,teq)%RHS(1,1,indx)   = dTdr_Top*sqrt(4.0D0*Pi)
                Endif

                If (fix_tvar_bottom) Then
                    equation_set(1,teq)%RHS(N_R,1,indx) = T_Bottom*sqrt(4.0D0*Pi)
                Endif


                If (fix_dtdr_bottom) Then
                    equation_set(1,teq)%RHS(N_R,1,indx) = dTdr_Bottom*sqrt(4.0D0*Pi)
                Endif

                If (magnetism) Then
                    ! No ell=0 equations for the B-field potential functions
                    equation_set(1,aeq)%RHS(:,:,indx)   = Zero
                    equation_set(1,ceq)%RHS(:,:,indx)   = Zero
                Endif

            Endif




            ! COMPRESSIBLE NOTE:   Leaving this for now.   We aren't using FE mode at the moment.
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
    End Subroutine Enforce_Boundary_Conditions

End Module Sphere_Linear_Terms
