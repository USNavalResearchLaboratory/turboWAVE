from setuptools import setup
from setuptools import find_packages

with open('README.md','r',encoding='utf-8') as fh:
	long_description = fh.read()

setup(name='twutils',
	version='5.0.0',
	license='PUBLIC-DOMAIN',
	author='Daniel Gordon',
	author_email='daniel.gordon@nrl.navy.mil',
	description='TurboWAVE utilities',
	long_description=long_description,
	long_description_content_type='text/markdown',
	url='https://github.com/USNavalResearchLaboratory/turboWAVE',
	classifiers=[
		'Programming Language :: Python :: 3',
	],
	package_dir={'':'src'},
	packages=find_packages(where='src'),
	python_requires='>=3.12',
	install_requires=[
		'numpy>=2.3',
		'scipy>=1.16',
		'matplotlib>=3.10',
		'jupyter>=1.1',
		'pillow>=11.3',
		'h5py>=3.14',
	],
	include_package_data=True,
	entry_points = { 'console_scripts': [
		'twinstall=twutils.command_line.twinstall:main',
		'twplot=twutils.command_line.twplot:main',
		'twtest=twutils.command_line.twtest:main',
		'twscan=twutils.command_line.twscan:main',
		'os2tw=twutils.command_line.os2tw:main',
		'wigner=twutils.command_line.wigner:main'] },
)
