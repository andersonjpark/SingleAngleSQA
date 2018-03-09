!***********************************************************************
!
!     units_module
!
!***********************************************************************

      module units_module

      type units_in_cgs
        real(8) :: e,pi,sinsqthetw     !numbers
        real(8) :: fm,km,MeV,s         !units
        real(8) :: c,ec,ga,G,GMs_c2km,gv,kB,mb,Ms,mu,Rs !cgs constants
        real(8) :: e0,h,GF,me,Q        !non-cgs constants
      end type units_in_cgs

      type(units_in_cgs), parameter :: units = units_in_cgs(
     &  2.718281828459045,          !e
     &  3.1415926535898e+00,        !pi
     &  0.2325,                     !sinsqthetw
     &  1.e-13,                     !fm
     &  1.e+05,                     !km
     &  1.602192e-06,               !MeV
     &  1.,                         !s
     &  2.997924562e+10,            !c
     &  4.80324214e-10,             !ec [sqrt(erg*cm)]
     &  1.23,                       !ga
     &  6.672387286039493e-08,      !G
     &  1.47664,                    !GMs_c2km
     &  1.,                         !gv
     &  1.38066244e-16,             !kB
     &  1.674e-24,                  !mb
     &  1.989e+33,                  !Ms
     &  1.661e-24,                  !mu
     &  6.9599e+10,                 !Rs
     &  8.9,                        !e0 [MeV/baryon]
     &  4.1356943e-21,              !h [MeV*s]
     &  8.957e-44,                  !GF [MeV*cm^3]
     &  0.511,                      !me [MeV]
     &  1.2935)                     !Q=neutron-proton mass [MeV]

      end module units_module

!***********************************************************************
