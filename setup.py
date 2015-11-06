from setuptools import setup, find_packages

setup(name='bamsurgeon',
	version='1.0',
	author='Adam Ewing',
	license='MIT',
	scripts=['bin/addindel.py', 'addsnv.py', 'addsv.py'],
	packages=find_packages(),
	)
