from setuptools import setup

def readme():
	with open('README.rst') as f:
		return f.read()

setup(name='twutils',
	version='3.3',
	description='TurboWAVE utilities',
	long_description=readme(),
	classifiers=[
		'Development Status :: 5 - Production/Stable',
		'Programming Language :: Python :: 3',
		],
	author='Daniel Gordon',
	author_email='daniel.gordon@nrl.navy.mil',
	url='http://www.nrl.navy.mil',
	packages=['twutils','twutils/QO'],
	install_requires=['scipy',],
	include_package_data=True,
	zip_safe=False)
