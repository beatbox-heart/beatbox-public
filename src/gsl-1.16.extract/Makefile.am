## Process this file with automake to produce Makefile.in

# AUTOMAKE_OPTIONS = readme-alpha

# SUBDIRS = gsl utils sys test err const complex cheb block vector matrix permutation combination multiset sort ieee-utils cblas blas linalg eigen specfunc dht qrng rng randist fft poly fit multifit statistics siman sum integration interpolation histogram ode-initval ode-initval2 roots multiroots min multimin monte ntuple diff deriv cdf wavelet bspline doc
SUBDIRS = gsl utils sys err complex block vector matrix cblas blas linalg eigen 


SUBLIBS = block/libgslblock.la blas/libgslblas.la complex/libgslcomplex.la eigen/libgsleigen.la err/libgslerr.la linalg/libgsllinalg.la matrix/libgslmatrix.la  sys/libgslsys.la vector/libgslvector.la
# SUBLIBS = block/libgslblock.la blas/libgslblas.la bspline/libgslbspline.la complex/libgslcomplex.la cheb/libgslcheb.la dht/libgsldht.la diff/libgsldiff.la deriv/libgslderiv.la eigen/libgsleigen.la err/libgslerr.la fft/libgslfft.la fit/libgslfit.la histogram/libgslhistogram.la ieee-utils/libgslieeeutils.la integration/libgslintegration.la interpolation/libgslinterpolation.la linalg/libgsllinalg.la matrix/libgslmatrix.la min/libgslmin.la monte/libgslmonte.la multifit/libgslmultifit.la multimin/libgslmultimin.la multiroots/libgslmultiroots.la ntuple/libgslntuple.la ode-initval/libgslodeiv.la ode-initval2/libgslodeiv2.la permutation/libgslpermutation.la combination/libgslcombination.la multiset/libgslmultiset.la poly/libgslpoly.la qrng/libgslqrng.la randist/libgslrandist.la rng/libgslrng.la roots/libgslroots.la siman/libgslsiman.la sort/libgslsort.la specfunc/libgslspecfunc.la statistics/libgslstatistics.la sum/libgslsum.la sys/libgslsys.la test/libgsltest.la utils/libutils.la vector/libgslvector.la cdf/libgslcdf.la wavelet/libgslwavelet.la

noinst_HEADERS = gsl_math.h gsl_pow_int.h gsl_nan.h gsl_machine.h gsl_mode.h gsl_precision.h gsl_types.h gsl_version.h gsl_minmax.h gsl_inline.h

noinst_SCRIPTS = gsl-config

# pkgconfigdir = $(libdir)/pkgconfig
# pkgconfig_DATA= gsl.pc

CLEANFILES = gsl.pc gsl-config

EXTRA_DIST = autogen.sh gsl-config.in gsl.pc.in configure.ac THANKS BUGS gsl.spec.in gsl.m4 test_gsl_histogram.sh pkgconfig.test

noinst_LTLIBRARIES = libgsl.la
libgsl_la_SOURCES = version.c
libgsl_la_LIBADD = $(SUBLIBS)
libgsl_la_LDFLAGS = -version-info $(GSL_LT_VERSION)
# noinst_HEADERS = gsl_math.h templates_on.h templates_off.h build.h
noinst_HEADERS += gsl_math.h templates_on.h templates_off.h build.h permutation/gsl_permutation.h

MINGW32_HOST = @MINGW32_HOST@
if MINGW32_HOST
libgsl_la_LIBADD += cblas/libgslcblas.la
libgsl_la_LDFLAGS += -no-undefined
endif

# m4datadir = $(datadir)/aclocal
# m4data_DATA = gsl.m4

# noinst_PROGRAMS = gsl-randist gsl-histogram

# gsl_randist_SOURCES = gsl-randist.c
# gsl_randist_LDADD = libgsl.la cblas/libgslcblas.la

# gsl_histogram_SOURCES = gsl-histogram.c
# gsl_histogram_LDADD = libgsl.la cblas/libgslcblas.la

# check_SCRIPTS = test_gsl_histogram.sh pkgconfig.test
# TESTS = test_gsl_histogram.sh pkgconfig.test

#bin_PROGRAMS = main dummy
#dummy_SOURCES = version.c
#dummy_LDADD = $(SUBLIBS)
#main_SOURCES = version.c env.c
#main_LDADD = libgsl.la

edit = $(SED) \
	-e 's|@prefix[@]|$(prefix)|g' \
	-e 's|@exec_prefix[@]|$(exec_prefix)|g' \
	-e 's|@libdir[@]|$(libdir)|g' \
	-e 's|@includedir[@]|$(includedir)|g' \
	-e 's|@GSL_CFLAGS[@]|$(GSL_CFLAGS)|g' \
	-e 's|@GSL_LIBM[@]|$(GSL_LIBM)|g' \
	-e 's|@GSL_LIBS[@]|$(GSL_LIBS)|g' \
	-e 's|@LIBS[@]|$(LIBS)|g' \
	-e 's|@VERSION[@]|$(VERSION)|g'

gsl-config gsl.pc: Makefile 
	@rm -f $@ $@.tmp
	@$(edit) '$(srcdir)/$@.in' >>$@.tmp
	@chmod a-w $@.tmp
	@mv $@.tmp $@
	@echo creating $@

gsl-config: $(srcdir)/gsl-config.in
gsl.pc: $(srcdir)/gsl.pc.in

