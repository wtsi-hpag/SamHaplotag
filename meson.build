project('SamHaplotag', 'cpp',
    license: 'MIT',
    meson_version: '>=0.57.1',
    version: '0.0.3'
)

flags = ['-Ofast']
if get_option('buildtype') == 'debug'
    flags = ['-O0', '-g', '-DDEBUG']
endif
flags += ['-DPV=' + meson.project_version()]

thread_dep = dependency('threads')
test('test SamHaplotag', executable('SamHaplotag', 'SamHaplotag.cpp', dependencies : thread_dep, install : true, cpp_args : flags), args : '--help')
test('test 10xSpoof', executable('10xSpoof', '10xSpoof.cpp', dependencies : thread_dep, install : true, cpp_args : flags), args : '--help')
test('test 16BaseBCGen', executable('16BaseBCGen', '16BaseBCGen.cpp', dependencies : thread_dep, install : true, cpp_args : flags), args : '--help')
