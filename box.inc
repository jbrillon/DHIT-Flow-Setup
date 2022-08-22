*                          File : box.inc

        integer      MaxK, MaxIn, MaxG

        parameter(   MaxK = 32,
     ,               MaxIn = 1001,
     ,               MaxG = 128 )


          common /BOX/
     ,                      Pi_const,
     ,                      K_max,
     ,                      VolBox,
     ,                      SpecE,
     ,                      Ein,
     ,                      Kin,
     ,                      EnK,
     ,                      CtK,
     ,                      DivNorm,
     ,                      AmpNorm,
     ,                      Div,
     ,                      Amp,
     ,                      Vel,
     ,                      Vel2,
     ,                      Xg,
     ,                      Yg,
     ,                      Zg,
     ,                      Xs,
     ,                      Ys,
     ,                      Zs,
     ,                      Ni1,Ni2,Nj1,Nj2,Nk1,Nk2,
     ,                      Ngrd,
     ,                      NumK,
     ,                      NumKout


          integer
     ,                      Ni1,Ni2,Nj1,Nj2,Nk1,Nk2,
     ,                      Ngrd,
     ,                      NumK,
     ,                      NumKout

          real *8
     ,                      Pi_const,
     ,                      K_max,
     ,                      VolBox,
     ,                      SpecE(0:MaxK),
     ,                      Ein(MaxIn),
     ,                      Kin(MaxIn),
     ,                      EnK(0:MaxK),
     ,                      CtK(0:MaxK),
     ,                      DivNorm(0:MaxK),
     ,                      AmpNorm(0:MaxK,3),
     ,                      Div(0:(2*MaxK+1)**3,3),
     ,                      Amp(0:(2*MaxK+1)**3,3),
     ,                      Vel (MaxG,MaxG,MaxG,3),
     ,                      Vel2(MaxG,MaxG,MaxG,3),
     ,                      Xg(MaxG,MaxG,MaxG),
     ,                      Yg(MaxG,MaxG,MaxG),
     ,                      Zg(MaxG,MaxG,MaxG),
     ,                      Xs(MaxG,MaxG,MaxG),
     ,                      Ys(MaxG,MaxG,MaxG),
     ,                      Zs(MaxG,MaxG,MaxG)

