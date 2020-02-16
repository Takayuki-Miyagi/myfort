obj/angular_momentum_couplings.o : src/angular_momentum_couplings.f90 obj/functions_from_c.o obj/general.o 
obj/functions_from_c.o : src/functions_from_c.f90 
obj/general.o : src/general.f90 
obj/iteration_methods.o : src/iteration_methods.f90 obj/linear_algebra.o 
obj/linear_algebra.o : src/linear_algebra.f90 obj/mat_vec_complex.o obj/mat_vec_double.o obj/mat_vec_single.o obj/matrix_complex.o obj/matrix_double.o obj/matrix_single.o obj/vector_complex.o obj/vector_double.o obj/vector_single.o obj/sngl_dble_cmplx.o obj/general.o 
obj/mat_vec_complex.o : src/mat_vec_complex.f90 obj/matrix_complex.o obj/vector_complex.o obj/general.o 
obj/mat_vec_double.o : src/mat_vec_double.f90 obj/matrix_double.o obj/vector_double.o obj/general.o 
obj/mat_vec_single.o : src/mat_vec_single.f90 obj/matrix_single.o obj/vector_single.o obj/general.o 
obj/matrix_complex.o : src/matrix_complex.f90 obj/vector_complex.o obj/general.o 
obj/matrix_double.o : src/matrix_double.f90 obj/vector_double.o obj/general.o 
obj/matrix_single.o : src/matrix_single.f90 obj/vector_single.o obj/general.o 
obj/myfort.o : src/myfort.f90 obj/iteration_methods.o obj/linear_algebra.o obj/wave_functions.o obj/profiler.o obj/physics_constants.o obj/general.o obj/store_couplings.o obj/angular_momentum_couplings.o obj/functions_from_c.o 
obj/physics_constants.o : src/physics_constants.f90 obj/general.o 
obj/profiler.o : src/profiler.f90 obj/general.o 
obj/sngl_dble_cmplx.o : src/sngl_dble_cmplx.f90 obj/matrix_complex.o obj/vector_complex.o obj/matrix_double.o obj/vector_double.o obj/matrix_single.o obj/vector_single.o 
obj/store_couplings.o : src/store_couplings.f90 obj/angular_momentum_couplings.o obj/functions_from_c.o obj/profiler.o obj/general.o 
obj/vector_complex.o : src/vector_complex.f90 obj/general.o 
obj/vector_double.o : src/vector_double.f90 obj/general.o 
obj/vector_single.o : src/vector_single.f90 obj/general.o 
obj/wave_functions.o : src/wave_functions.f90 obj/general.o obj/physics_constants.o obj/functions_from_c.o 
