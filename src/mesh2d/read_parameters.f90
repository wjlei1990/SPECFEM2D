subroutine read_parameters(nex,nez,model_x1,model_x2,model_z1,model_z2,nproc_x,nproc_z,&
        simul_type,debug)
   integer ::nex,nez,nex_one_proc,nez_one_proc
   double precision ::model_x1,model_x2,model_z1,model_z2
   integer ::nproc_x,nproc_z
   integer ::simul_type
   logical ::debug

   nex=20
   nez=20
   nproc_x=5
   nproc_z=5
   nex_one_proc=nex/nproc_x
   nez_one_proc=nez/nproc_z
   model_x1=0
   model_x2=10000
   model_z1=0
   model_z2=10000
   simul_type=2
   debug=.false.

end subroutine
