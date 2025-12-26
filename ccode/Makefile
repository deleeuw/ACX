
SSRC = smacofIsotone.c smacofMPInverseV.c smacofSort.c smacofSSEngine.c smacofSSMajorize.c smacofSSMonotone.c

%.o: %.c smacofSS.h
	clang -c $@

shlib: smacofSS.h $(SSRC)
	R CMD SHLIB -o smacofSS.so $(SSRC)

clean:
	rm -f *.o

pristine:
	rm -f *.o *.so

