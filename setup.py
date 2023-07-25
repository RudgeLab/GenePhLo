from setuptools import setup, find_packages


with open('README.md', 'r') as ld:
    long_description = ld.read()

setup(name='genephlo',
    version='0.0.0',
    author='Tim Rudge',
    author_email='timrudge@gmail.com',
    description='GenePhLo: Genetic Phase-based Logic',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/RudgeLab/GenePhLo',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        ],
    install_requires=[        
	'loica',
        'numpy', 
	'scipy',
    	'pandas',
        'matplotlib'
        ],
    setup_requires=[
        'loica',
	'numpy',
	    'scipy',
	    'pandas',
            'matplotlib'
        ],
    packages=find_packages(),
    python_requires='>=3'
)
