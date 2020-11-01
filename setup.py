import setuptools

with open('README.md','r') as fh:
	long_description = fh.read()
	
setuptools.setup(
	name='POSMM',
	packages=setuptools.find_packages(),
	entry_points={
		'console_scripts':[
			'POSMM = POSMM.POSMM:main'
		]
	},
	package_data={'POSMM':['posmmbin/*','posmmsource/*cpp','models/*json']},
	version='1.0',
	author='David Burks',
	author_email='davidburks@my.unt.edu',
	description='Python-Optimized Single-Order Markov Model Classifier for Metagenomic Reads',
	long_description=long_description,
	url='https://github.com/djburks/POSMM',
	install_requires=[
		'sklearn_json',
		'sklearn',
		'ncbitax2lin',
	],
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: MIT License",
		"Operating System :: GNU/Linux",
	],
	python_requires='>3.7',
)
	 
