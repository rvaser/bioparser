project(
  'bioparser',
  ['cpp'],
  version : '3.1.0',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++11',
    'b_ndebug=if-release'],
  license : 'MIT',
  meson_version : '>=0.60.0'
)

################
# Dependencies #
################

# zlib
bioparser_zlib_dep = dependency('zlib', include_type : 'system', fallback : ['zlib', 'zlib_dep'])

bioparser_lib_deps = [
  bioparser_zlib_dep,
]

###########
# Headers #
###########

subdir('include')

#########
# Tests #
#########

if (not meson.is_subproject()) and get_option('tests')
  # biosoup
  bioparser_biosoup_dep = dependency('biosoup', version : '>= 0.11.0', fallback : ['biosoup', 'biosoup_dep'])

  # gtest
  bioparser_gtest_dep = dependency('gtest', version : '>= 1.10.0', main : true, fallback : ['gtest', 'gtest_main_dep'])

  subdir('test')
endif

###################
# Dependency info #
###################

if (not meson.is_subproject())
  import('pkgconfig').generate(
    name : 'bioparser',
    version : meson.project_version(),
    filebase : 'bioparser',
    description : 'C++ header only parsing library for bioinformatics formats (FASTA/Q, MHAP/PAF/SAM), with support for zlib compressed files.')
endif

bioparser_dep = declare_dependency(
  include_directories : bioparser_include_directories,
  dependencies : bioparser_lib_deps,
  version : meson.project_version()
)
