    Character(len=128) :: build_git_commit = "---git ha"// &
        "sh---e95c9d10165d8454cdfdd7f2e448930edba09835 +++ Tue A"// &
        "ug 20 10:37:42 2019 -0600"//"---"//Char(0)

    Character(len=128) :: build_git_url = "---git url--"// &
        "-git@github.com:geodynamics/Rayleigh.git"//"---"//Char(0)

    Character(len=128) :: build_git_branch = "---git br"// &
        "anch---master"//"---"//Char(0)

    Character(len=128) :: build_date = "---build date--"// &
        "-Mon Sep  9 09:16:42 MDT 2019"//"---"//Char(0)

    Character(len=256) :: build_machine = "---build mac"// &
        "hine---Linux sanitas 4.15.0-60-generic #67-Ubuntu SMP T"// &
        "hu Aug 22 16:55:30 UTC 2019 x86_64 x86_64 x86_64 GNU/Li"// &
        "nux"//"---"//Char(0)

    Character(len=512) :: build_dir = "---build dir---/"// &
        "home/feathern/devel/forks/Rayleigh"//"---"//Char(0)

    Character(len=512) :: build_custom_dir = "---custom"// &
        " dir---/home/feathern/devel/racomp/patch_nick2"//"---"//Char(0)

    Character(len=64 ) :: build_version = "---build ver"// &
        "sion---dbg"//"---"//Char(0)

    Character(len=512) :: build_fflags = "---build ffla"// &
        "gs----O0 -g -fbounds-check -ffpe-trap=invalid,zero,over"// &
        "flow -fbacktrace -ffixed-line-length-132 -DGNU_COMPILER"//"---"//Char(0)

    Character(len=512) :: build_lib = "---build lib----"// &
        "L/custom/software/miniconda3/envs/radev/lib -Wl,--no-as"// &
        "-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpth"// &
        "read -lm -ldl -lstdc++"//"---"//Char(0)

    Character(len=512) :: FC_version = "---FC version--"// &
        "-GNU Fortran (crosstool-NG 1.23.0.449-a04d0) 7.3.0"//"---"//Char(0)

    Character(len=512) :: FC_location = "---FC location"// &
        "---/custom/software/miniconda3/envs/radev/bin/mpifort"//"---"//Char(0)

