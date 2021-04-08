project('SamHaplotag', 'cpp',
    license: 'MIT',
    meson_version: '>=0.57.1',
    version: '0.0.1'
)

thread_dep = dependency('threads')
executable('SamHaplotag', 'SamHaplotag.cpp', dependencies : thread_dep, install : true)
executable('10xSpoof', '10xSpoof.cpp', dependencies : thread_dep, install : true)

Test = find_program('Test/run.sh')
test('Test', Test)
