bioparser_test_config = configuration_data()
bioparser_test_config.set('BIOPARSER_TEST_DATA',  meson.project_source_root() + '/test/data/')

###########
# Sources #
###########

bioparser_test_config_headers = [
  configure_file(
    input : 'bioparser_test_config.h.in',
    output : 'bioparser_test_config.h',
    configuration : bioparser_test_config)
]

bioparser_test_sources = files([
  'fasta_parser_test.cpp',
  'fastq_parser_test.cpp',
  'mhap_parser_test.cpp',
  'paf_parser_test.cpp',
  'parser_test.cpp',
  'sam_parser_test.cpp',
]) + bioparser_test_config_headers

bioparser_test = executable(
  'bioparser_test',
  bioparser_test_sources,
  dependencies : [bioparser_lib_deps, bioparser_biosoup_dep, bioparser_gtest_dep],
  include_directories : bioparser_include_directories,
  install : false)

#########
# Tests #
#########

test('bioparser gtest unit tests', bioparser_test)
