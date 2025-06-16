#
# ParalinearDesign
#
# Copyright 2025 Seoul National University
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# “Software”), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
# NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

import os
import sys
import shutil
from setuptools import setup
from setuptools.command.build_py import build_py
from distutils import ccompiler
import distutils.sysconfig


PLD_BINARY_NAME = 'paralinear_design'


class BuildPyWithParalinearDesignBinary(build_py):

    def run(self):
        cpp_flags = ['-std=c++11', '-DSPECIAL_HP', '-DPLD_CAPABILITY=2']
        source_file = 'src/linear_design.cpp'
        build_cmd = self.get_finalized_command('build')
        output_dir = 'pldesign'
        output_binary_name = PLD_BINARY_NAME
        output_binary_path = os.path.join(output_dir, output_binary_name)

        library_path = os.path.join('src', 'Utils', 'libraries', 'LinearDesign_linux64.so')
        include_dirs = ['src', 'src/Utils']

        os.makedirs(output_dir, exist_ok=True)
        build_temp_dir = build_cmd.build_temp

        os.makedirs(build_temp_dir, exist_ok=True)
        compiler_type = build_cmd.compiler
        compiler = ccompiler.new_compiler(compiler=compiler_type)
        distutils.sysconfig.customize_compiler(compiler)

        if not hasattr(compiler, 'compiler_cxx'):
             compiler.compiler_cxx = [compiler.compiler[0]]

        try:
            objects = compiler.compile(
                [source_file],
                output_dir=build_temp_dir,
                include_dirs=include_dirs,
                extra_preargs=cpp_flags
            )
        except Exception as e:
            print(f'Compilation failed: {e}', file=sys.stderr)
            sys.exit(1)

        lib_name = os.path.basename(library_path)
        if lib_name.startswith('lib'):
            lib_name = lib_name[3:]
        if lib_name.endswith('.so'):
            lib_name = lib_name[:-3]
        elif lib_name.endswith('.dylib'):
            lib_name = lib_name[:-6]

        shutil.copy(library_path,
                    os.path.join(output_dir,
                                 'lib' + os.path.basename(library_path)))

        try:
            compiler.link_executable(
                objects,
                output_binary_name,
                output_dir=output_dir,
                libraries=[lib_name],
                library_dirs=[output_dir],
                extra_postargs=['-Wl,-rpath=$ORIGIN'],
                debug=build_cmd.debug,
                target_lang='c++'
            )
        except Exception as e:
            print(f'Linking failed: {e}', file=sys.stderr)
            sys.exit(1)

        try:
            os.chmod(output_binary_path, os.stat(output_binary_path).st_mode | 0o111)
            print(f'Successfully compiled and linked {output_binary_path}')
        except OSError as e:
            print(f'Error setting executable permission: {e}', file=sys.stderr)

        # Copy the data files to the package directory
        data_files = [
            'coding_wheel.txt',
            'codon_usage_freq_table_human.csv',
            'codon_usage_freq_table_yeast.csv',
        ]
        for data_file in data_files:
            shutil.copy(data_file, os.path.join(output_dir, data_file))

        super().run()

setup(
    name='pldesign',
    version='0.1',
    description='Codon Optimizer for mRNA Vaccine Design',
    author='Hyeshik Chang',
    author_email='hyeshik@snu.ac.kr',
    url='https://github.com/ChangLabSNU/ParalinearDesign',
    download_url='https://github.com/ChangLabSNU/ParalinearDesign/releases',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    keywords=[
        'mRNA vaccine',
        'messenger RNA',
        'codon optimization'
    ],
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=['pldesign'],
    package_data={
        'pldesign': [
            'coding_wheel.txt',
            'codon_usage_freq_table_human.csv',
            'codon_usage_freq_table_yeast.csv',
            PLD_BINARY_NAME,
            'libLinearDesign_linux64.so',
        ],
    },
    entry_points={
        'console_scripts': [
            'pldesign = pldesign.main:main',
        ],
    },
    install_requires=[
        'numpy >= 2.0',
        'pandas >= 2.0',
        'tqdm >= 4.0',
        'biopython >= 1.5',
        'viennarna >= 2.6.2'
    ],
    cmdclass={
        'build_py': BuildPyWithParalinearDesignBinary
    },
    include_package_data=True,
    zip_safe=False,
)
