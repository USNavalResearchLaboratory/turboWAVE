from setuptools import setup

def readme():
	with open('README.rst') as f:
		return f.read()

setup(name='twutils',
	version='4.0.1',
	description='TurboWAVE utilities',
	long_description=readme(),
	classifiers=[
		'Development Status :: 5 - Production/Stable',
		'Programming Language :: Python :: 3',
		],
	author='Daniel Gordon',
	author_email='daniel.gordon@nrl.navy.mil',
	url='https://github.com/USNavalResearchLaboratory/turboWAVE',
	packages=['twutils','twutils/QO'],
	install_requires=['scipy',],
	include_package_data=True,
	zip_safe=False)
