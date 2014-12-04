clean:
	rm -rf *.so *~
simple:
	f2py -m add add.f

comp:
	f2py -c -m add add.f

better:
	f2py -h add.pyf -m add add.f
	f2py -c -m add add.f add.pyf
