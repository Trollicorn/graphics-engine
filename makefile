run: main.py scripts/test2.mdl
	python3 main.py scripts/run.mdl

art: main.py scripts/art.mdl
	python3 main.py scripts/art.mdl

flat: main.py scripts/test4.mdl
	python3 main.py scripts/test4.mdl

sphere: main.py scripts/test3.mdl
	python3 main.py scripts/test3.mdl

sflat: main.py scripts/test5.mdl
	python3 main.py scripts/test5.mdl

obj: main.py scripts/test6.mdl
	python3 main.py scripts/test6.mdl

clean:
	rm */*pyc */*out */parsetab.py
	rm -rf */__pycache__

clear:
	rm */*pyc *out */parsetab.py */*ppm */*png
	rm -rf */__pycache__
