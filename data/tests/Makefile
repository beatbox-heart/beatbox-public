#COMMENT=@echo making "'$@'" stem "'$*'" out of "'$+'" due to "'$?'"
COMMENT=\#@ignore
BBXS=Beatbox_SEQ
BBXP=Beatbox
MPIOPT:=${shell mpirun --version | perl -ane '$$.==1 and $$v=pop @F and $$f=substr($$v,0,1) and $$f>=3 and print("--oversubscribe");'}
MPIRUN:=export TMPDIR=/tmp; mpirun $(MPIOPT) -np

TEST=echo -n test $* :
OK=echo ' ok'
FAIL=echo ' failed'; ./admitfailure

NEGATIVE=nodev depth nexist notype reserved declared diff0d nomem nomem-4
POSITIVE=fhn0 fhn1 fhn2 fhn3 fhn3side fhn3side-16 bkl3cont-0 bkl3cont-4 br0 reduce reduce-4 hh52mrl
ALL=$(NEGATIVE) $(POSITIVE)
all: closeenough ignore admitfailure ${foreach t,$(ALL),$t.test}
ignore: ignore.c
	cc $< -o $@

.PHONY: %.test

##########################
# Negative tests: scripts leading to exceptions; check for error messages in the log files.
%.test: %.log %.errmsg
	$(COMMENT)
	@$(TEST)
	@if grep "`cat $*.errmsg`" $< > /dev/null; then $(OK); else $(FAIL); fi
	@rm $*.log
%.log: %.bbs
	@-$(BBXS) $< -log $@ > /dev/null 2>/dev/null

%-4.log: %.bbs
	@-$(MPIRUN) 4 $(BBXP) $< > $@ 2>/dev/null

##########################
# Positive tests: output files compared to what they should be: either exact or approximate. 
%.test: %.out %.exact
	$(COMMENT)
	@$(TEST)
	@if cmp $+; then $(OK); else $(FAIL); fi
	@rm $*.out
%.test: closeenough %.out %.approx
	$(COMMENT)
	@$(TEST)
	@if ./$+; then $(OK); else $(FAIL); fi
	@rm $*.out
# Outputs of parallel jobs may come out in arbitrary order
%.test: closeenough %.sorted %.sortedapprox
	$(COMMENT)
	@$(TEST)
	@if ./$+; then $(OK); else $(FAIL); fi
	@rm $*.sorted
	@rm $*.out
%.sorted: %.out
	$(COMMENT) '%.sorted'
	@sort $< > $@
# .PRECIOUS: fhn0.out
%.out: %.bbs
	$(COMMENT) '%.out'
	@-rm -rf $@
	@$(BBXS) $< $@ -log /dev/null > /dev/null
%.fullout: %.bbs
	$(COMMENT) '%.fullout'
	$(BBXS) $< $@ -log $@.log -verbose | tee $@.stdout
%-4.out: %.bbs
	$(COMMENT)
	@-rm -rf $@
	@$(MPIRUN) 4 $(BBXP) $< $@ -log /dev/null > /dev/null
%-4.fullout: %.bbs
	$(COMMENT)
	@-rm -rf $@
	@$(MPIRUN) 4 $(BBXP) $+ $@ -log $@.log -verbose > $@.stdout
%-16.out: %.bbs
	$(COMMENT)
	@-rm -rf $@
	@$(MPIRUN) 16 $(BBXP) $+ $@ -decomp 4x4x1 -log /dev/null > /dev/null
%-16.fullout: %.bbs
	$(COMMENT)
	@-rm -rf $@
	@$(MPIRUN) 16 $(BBXP) $+ $@ -decomp 4x4x1 -log /dev/null -verbose > $@.stdout
fhn1.out: fhn0.out fhn.par
fhn2.out: fhn1.out fhn.par
fhn3.out: fhn1.out fhn.par

bkl3.half bkl3.full: bkl3.bbs
	$(COMMENT)
	@$(BBXS) $< bkl3.half bkl3.full  -log /dev/null > /dev/null
.PRECIOUS: bkl3cont-%.exact
bkl3cont-%.exact: bkl3.full
	$(COMMENT)
	@cp -p $< $@
bkl3cont-0.out: bkl3cont.bbs bkl3.half
	@$(BBXS) $+ $@ -log /dev/null > /dev/null
bkl3cont-4.out: bkl3cont.bbs bkl3.half
	@$(MPIRUN) 4 $(BBXP) $+ $@ -log /dev/null > /dev/null

##########################
# Auxiliary utilities
closeenough: closeenough.c
	$(COMMENT)
	cc $< -o $@
admitfailure: admitfailure.c
	$(COMMENT)
	cc $< -o $@
