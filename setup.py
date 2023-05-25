from setuptools import setup, find_packages

setup(
    name='jump_detection',
    # version='0.1.0',
    # url='https://github.com/yourname/my_project',
    # author='Author Name',
    # author_email='author@gmail.com',
    # description='Description of my package',
    packages=find_packages(),    
    install_requires=['numpy >= 1.18.1', 'pandas >= 1.0.1', 'scipy >= 1.4.1'],
)
