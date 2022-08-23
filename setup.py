from setuptools import setup, find_packages, Extension

setup(
    name="mykmeanspp",
    version="0.0.1",
    author="Roei Ben Zion",
    author_email='roei.benzion@gmail.com',
    description='mykmeanssp',
    install_requires=['invoke'],
    py_modules=['mykmeanssp'],
    packages=find_packages(),  # find_packages(where='.', exclude=())
                               #    Return a list of all Python packages found within directory 'where'
    ext_modules=[
        Extension(
            'mykmeanssp', 
            ['kmeans.c'],
            ),
            ]
    )
