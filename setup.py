from setuptools import setup
import io
import os

import codonOptim

here = os.path.abspath(os.path.dirname(__file__))

def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.md')

setup(
    name='codonOptim',
    version=codonOptim.__version__,
		url='https://github.com/haydnKing/codon-optim',
    license='Apache Software License',
    author='Haydn King',
    install_requires=['Bio>=1.65',
											'numpy>=1.11.0',
                    ],
    cmdclass={'test': PyTest},
    author_email='hjking734@gmail.com',
    description='Smart Codon Optimisation',
    long_description=long_description,
    packages=['codonOptim'],
    include_package_data=True,
    platforms='any',
)
