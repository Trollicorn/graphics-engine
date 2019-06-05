run: main.py scripts/test.mdl
	python3 main.py scripts/test.mdl

clean:
	rm */*pyc */*out */parsetab.py
	rm -rf */__pycache__

clear:
	rm */*pyc *out */parsetab.py */*ppm */*png
	rm -rf */__pycache__
