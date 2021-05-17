from setuptools import setup

install_requires = [
    'matplotlib'
    ]

setup(
    name='TASplot',
    version=0.1,
    description="""
        tasplot is a Python module that adds the names of volcanic rock
        types to a plot of total alkali versus silica data""",
    license='MIT',
    author='John A. Stevenson',
    author_email='@volcan01010',
    install_requires=install_requires,
    packages=['tasplot'],
    classifiers=[
        'Development Status :: 3 - Beta',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python'
    ]
)
