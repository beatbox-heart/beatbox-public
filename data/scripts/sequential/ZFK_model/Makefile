# Testing various solvers using analytical solution of the ZFK 
# (Zel'dovich - Frank-Kamenetsky model, a.k.a. Nagumo equation

# Make a square geometry NxN for a given N
square_%.bbg: square.pl
	@$< $* > $@

zfkmon2.dir: zfkmon2.bbs square_20.bbg
	$(BBXS) $<

zfkbid2.dir: zfkbid2.bbs square_20.bbg
	$(BBXS) $<
