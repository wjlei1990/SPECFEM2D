module simulation_type

  use wave2d_variables

contains

  !---------------------------------------------
  subroutine set_simulation_flag()
    
    if(SIMUL_TYPE==1) then
      !dynamic relaxation
      SOURCE_FLAG=.false.
      FORCE_FLAG=.true.
      DAMP_FLAG=.true.
      ABSORB_FLAG=.false.
      S_ALPHA=0.5
      S_BETA=0.
    elseif(SIMUL_TYPE==2) then
      !wave simulation
      SOURCE_FLAG=.true.
      FORCE_FLAG=.false.
      DAMP_FLAG=.false.
      ABSORB_FLAG=.true.
      if(IM_TRUE)then
        S_ALPHA=0.5
        S_BETA=0.25
      else
        S_ALPHA=0.5
        S_BETA=0.
      endif
    endif

  end subroutine set_simulation_flag

end module simulation_type

