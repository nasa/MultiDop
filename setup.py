from setuptools import setup
import subprocess

setup(
    name='MultiDop',
    version='0.3',
    author='Timothy Lang',
    author_email='timothy.j.lang@nasa.gov',
    packages=['multidop',],
    license='LICENSE.md',
    description='Multiple-Doppler Radar Analysis Toolkit (MultiDop)',
    long_description=open('description.txt').read(),
    install_requires=['arm_pyart'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console"
        ],
)

subprocess.call('cd src; make', shell=True)
