cd Meschach
del *.obj meschach.lib
cl /c /Ox /Op /DNDEBUG copy.c err.c matrixio.c memory.c vecop.c matop.c pxop.c 	submat.c init.c otherio.c machine.c matlab.c ivecop.c version.c meminfo.c memstat.c lufactor.c bkpfacto.c chfactor.c qrfactor.c solve.c hsehldr.c givens.c update.c norm.c hessen.c symmeig.c schur.c svd.c fft.c mfunc.c bdfactor.c
link /lib /out:meschach.lib copy.obj err.obj matrixio.obj memory.obj vecop.obj matop.obj pxop.obj 	submat.obj init.obj otherio.obj machine.obj matlab.obj ivecop.obj version.obj meminfo.obj memstat.obj lufactor.obj bkpfacto.obj chfactor.obj qrfactor.obj solve.obj hsehldr.obj givens.obj update.obj norm.obj hessen.obj symmeig.obj schur.obj svd.obj fft.obj mfunc.obj bdfactor.obj
cd ..

cd brent
del *.obj brent.lib
cl /c /Ox /Op /DNDEBUG fminbr.c
link /lib /out:brent.lib fminbr.obj
cd ..

cd ..\libdict-0.2.1
del *.obj dict.lib
cl /c /Ox /Op /DNDEBUG dict.c hashtable.c hb_tree.c pr_tree.c rb_tree.c sp_tree.c tr_tree.c wb_tree.c
link /lib /out:dict.lib dict.obj hashtable.obj hb_tree.obj pr_tree.obj rb_tree.obj sp_tree.obj tr_tree.obj wb_tree.obj
cd ..\src

del *.obj slr.exe
cl /Feslr /Ox /Op /DWINDOWS /DNDEBUG statistics.c gamma.c bases.c codonmodel.c data.c gencode.c like.c matrix.c model.c optimize.c options.c rng.c slr.c linemin.c spinner.c tree.c tree_data.c utility.c math_win.c  mystring.c Meschach\meschach.lib brent\brent.lib ..\libdict-0.2.1\dict.lib