from setuptools import setup, find_packages

version = '1.0'

setup(name='mifish',
      version=version,
      description="the command line version of MiFish pipeline. It can also be used with any other eDNA meta-barcoding primers",
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
          'Operating System :: Unix',
          'Programming Language :: Python :: 3.9',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      keywords='mifish',
      author='Tao Zhu',
      author_email='zhutao@edu.k.u-tokyo.ac.jp',
      url='http://mitofish.aori.u-tokyo.ac.jp/mifish',
      license='GPLv3+',
      packages=find_packages(),
      include_package_data=True,
      python_requires='>=3.8',
      entry_points={
          'console_scripts': [
            'mifish = mifish.cmd.mifish:main',
          ]
      },
)