# Scalability test.

# Shorthand for both version of Beatbox
BBXP=Beatbox
BBXS=Beatbox_SEQ

scroll-%.out: scroll.bbs bkl1.rec
	( time mpirun -np $* $(BBXP) $+ $* ) >$@ 2>&1


# Sequential run for comparison
scroll.out: scroll.bbs bkl1.rec
	( time $(BBXS) $+ 0 ) >$@ 2>&1

ALLPAR=${foreach n,1 2 4 8 16,scroll-$n.out}
all: scroll.out $(ALLPAR)

report.dat: $(ALLPAR)
	grep took $(ALLPAR) | tee $ | perl -pe 's/scroll-//;s/.out:This run took/\t/;s/seconds.//' | tee $@

# Visualization of the result. 
# Uses ezview v. 1.1, see http://empslocal.ex.ac.uk/people/staff/vnb262/software/EZView/
.PHONY: view-%
view-%: view.pl scroll-%.dir
	$+
