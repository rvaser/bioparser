###########
# Headers #
###########

bioparser_include_directories = include_directories(['.'])

if meson.is_subproject()
  subdir_done()
endif

install_headers(
  files([
    'bioparser/fasta_parser.hpp',
    'bioparser/fastq_parser.hpp',
    'bioparser/mhap_parser.hpp',
    'bioparser/paf_parser.hpp',
    'bioparser/parser.hpp',
    'bioparser/sam_parser.hpp',
  ]),
  subdir : 'bioparser')
