from setuptools import setup

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='ksp-astrogator',
    version='0.0.1',
    packages=[''],
    url='https://github.com/msparapa/ksp-astrogator',
    license='',
    author='Michael Sparapany',
    author_email='',
    description='',
    install_requires=requirements
)
