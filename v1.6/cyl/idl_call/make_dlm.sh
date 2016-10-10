ftnchek -nocheck -quiet -makedcl galaxy_idl.for
\rm galaxy_idl_load.c
dlmform galaxy_idl.for
\rm galaxy_idl.f
