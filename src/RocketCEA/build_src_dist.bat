

cd .\docs
D:\python37_64\Scripts\sphinx-build.exe -b html -d _build/doctrees  . _build/html
cd ..

xcopy /S .\docs\_build\html .\rocketcea\sphinx_html /y

python setup.py sdist
