project('climate model', 'c', 'fortran',
	version: '1.0.0')

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)
gf_dep = cc.find_library('gfortran', required : true)
gsl_dep = cc.find_library('gsl', required : true)
gslcblas_dep = cc.find_library('gslcblas', required : true)
gsl_deps = [gsl_dep, gslcblas_dep]

rrtm_deps = [gf_dep]
foreach libname: ['libfpda_rrtm_sw', 'libfpda_rrtm_lw']
    rrtm_deps += [cc.find_library(libname, dirs: [meson.source_root() + '/' + 'rrtm/install/lib/'])]
endforeach

# Less precise but faster calculations
c_flags = ['-ffast-math']

rrtm_incdir = include_directories('rrtm/install/include/')
gnuplot_i_incdir = include_directories('gnuplot_i/src/')
ascii_incdir = include_directories('ascii/')
eddington_incdir = include_directories('eddington/')

ascii_lib = shared_library('ascii', 'ascii/ascii.c', version: '1.0.0', soversion: '0', dependencies: m_dep)
eddington_lib = shared_library('eddington', 'eddington/eddington_v2.c', version: '2.0.0', soversion: '0', dependencies: m_dep)
gnuplot_i_lib = shared_library('gnuplot_i', 'gnuplot_i/src/gnuplot_i.c', version: '1.0.0', soversion: '0')
rt_lib = shared_library('thermal_radiation', 'thermal_radiation.c', link_with: [eddington_lib], include_directories: [eddington_incdir], version: '1.0.0', soversion: '0', dependencies: m_dep)

executable('model', 'model.c', dependencies : [m_dep] + rrtm_deps + gsl_deps, link_with: [rt_lib, ascii_lib, gnuplot_i_lib], include_directories: [ascii_incdir, gnuplot_i_incdir, gnuplot_i_incdir, rrtm_incdir], c_args: c_flags)
executable('opacity-analysis', 'opacity-analysis.c', dependencies : m_dep, link_with: rt_lib)

e = executable('thermal_radiation_test', 'thermal_radiation_test.c', dependencies: m_dep, link_with: rt_lib)
test('Thermal Radiation Test', e)
