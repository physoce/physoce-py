from setuptools import setup

# read the contents of README file
# see https://packaging.python.org/guides/making-a-pypi-friendly-readme/
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='physoce',
      version='0.0.2',
      url='https://github.com/physoce/physoce-py',
      description='Python tools for Physical Oceanography',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Tom Connolly',
      author_email='tconnolly@mlml.calstate.edu',
      packages=['physoce'])
