from setuptools import setup, find_packages

setup(name='bamsurgeon',
	version='1.0',
	author='Adam Ewing',
	license='MIT',
	scripts=['bin/*']
	packages=find_packages(),
	)
