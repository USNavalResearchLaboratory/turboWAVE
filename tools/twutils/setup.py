from setuptools import setup
from setuptools import find_packages

def readme():
	with open('README.rst') as f:
		return f.read()

setup(name='twutils',
	version='4.2.1',
	description='TurboWAVE utilities',
	long_description=readme(),
	classifiers=[
		'Development Status :: 5 - Production/Stable',
		'Programming Language :: Python :: 3',
		],
	python_requires='>=3.8',
	author='Daniel Gordon',
	author_email='daniel.gordon@nrl.navy.mil',
	url='https://github.com/USNavalResearchLaboratory/turboWAVE',
	packages=find_packages(),
	install_requires=['scipy>=1.5','matplotlib>=3.3','jupyter>=1','ipympl>=0.5','pillow>=8','h5py'],
	include_package_data=True,
	entry_points = { 'console_scripts': [
		'twplot=twutils.command_line.twplot:main',
		'twtest=twutils.command_line.twtest:main',
		'os2tw=twutils.command_line.os2tw:main',
		'wigner=twutils.command_line.wigner:main'] },
	zip_safe=False)
