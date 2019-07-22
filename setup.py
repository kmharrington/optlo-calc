import os
from setuptools import setup, Command, find_packages

def readme():
    with open('README.rst') as f:
        return f.read()

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf build dist *.pyc *.tgz *.egg-info')


setup(name = 'optlo_calc',
      version = '0.1',
      description = '''Code for Calculating Optical Loading or Transmission 
                       through microwave optical chains''',
      long_description = readme(),
      author = 'Katie Harrington',
      author_email = 'katie.megan.harrington@gmail.com',
      license = 'MIT',
      packages = ['optlo_calc', 'optlo_calc.optics'],
      package_dir = {'optlo_calc':'python'},
      cmdclass={'clean':CleanCommand,},
      #scripts = []
      )