obj/angular_momentum_couplings.o : src/angular_momentum_couplings.f90 obj/functions_from_c.o 
obj/functions_from_c.o : src/functions_from_c.f90 
obj/general.o : src/general.f90 
obj/iteration_methods.o : src/iteration_methods.f90 obj/linear_algebra.o 
obj/linear_algebra.o : src/linear_algebra.f90 obj/matrix_definitions.o obj/vector_definitions.o 
obj/matrix_definitions.o : src/matrix_definitions.f90 obj/vector_definitions.o 
obj/myfort.o : src/myfort.f90 obj/renormalization.o obj/iteration_methods.o obj/linear_algebra.o obj/wave_functions.o obj/profiler.o obj/physics_constants.o obj/general.o obj/store_couplings.o obj/angular_momentum_couplings.o obj/functions_from_c.o 
obj/physics_constants.o : src/physics_constants.f90 
obj/profiler.o : src/profiler.f90 obj/general.o 
obj/renormalization.o : src/renormalization.f90 obj/profiler.o obj/linear_algebra.o 
obj/store_couplings.o : src/store_couplings.f90 obj/angular_momentum_couplings.o obj/functions_from_c.o obj/profiler.o 
obj/vector_definitions.o : src/vector_definitions.f90 obj/general.o 
obj/wave_functions.o : src/wave_functions.f90 obj/physics_constants.o obj/functions_from_c.o 
obj/dvode_f90_m.o : src/dvode/dvode_f90_m.f90
obj/renormalization.o : src/renormalization.f90 obj/dvode_f90_m.o obj/linear_algebra.o obj/profiler.o
