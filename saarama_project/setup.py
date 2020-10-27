from setuptools import setup, find_packages

setup(
      # mandatory
      name='saarama',
      # mandatory
      version='0.1',
      # mandatory
      author_email='joshua.bopp@stud-mail.uni-wuerzburg.de',
      packages=['saarama'],
      package_data={},
      install_requires=['matplotlib', 'MDAnalysis', 'click', 'Biopython', 'seaborn', 'scipy'],
      entry_points={
        'console_scripts': ['saarama = saarama.cli:start']
      }
)
